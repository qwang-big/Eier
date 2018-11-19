importBW <- function(meta, granges, threads=1) {
  files <- if (class(meta)=="matrix"|class(meta)=="data.frame")
    meta[,'file']
  else
    meta
  files <- as.character(files)
  gr <- if (class(granges)!="data.frame")
    as.data.frame(granges)
  else
    granges
  if (threads>1) {
    matrix(unlist(mclapply(files, function(file){
    cat(paste0("reading ", file, " ...\n"))
    mclapply(gr, 1, function(r) {
      readBW(file, r[1], as.integer(r[2]), as.integer(r[3]))
    })}, mc.cores = threads)), ncol=length(meta[,'file']), byrow=FALSE)
  } else {
    matrix(unlist(lapply(files, function(file){
    cat(paste0("reading ", file, " ...\n"))
    apply(gr, 1, function(r) {
      readBW(file, r[1], as.integer(r[2]), as.integer(r[3]))
    })})), ncol=length(meta[,'file']), byrow=FALSE)
  }
}

readBW <- function(bwfile, chrom, start, end) {
  read_bigwig_impl(bwfile, chrom, start, end)
}

readBED <- function(file){
  if (class(file)!="data.frame")
    bed=read.table(file,stringsAsFactors = FALSE)
  else
    bed=file
  if (ncol(bed)==4)
    GRanges(Rle(bed[,1]), IRanges(as.integer(bed[,2]), as.integer(bed[,3])), id=bed[,4])
  else
    GRanges(Rle(bed[,1]), IRanges(as.integer(bed[,2]), as.integer(bed[,3])), id=paste0(bed[,1],'_',bed[,2],'_',bed[,3]))
}

filterPeak <- function(files, ref, group=NULL, filter = c("both", "enhancer", "promoter")){
  filter = match.arg(filter)
  i <- c()
  if (class(ref)=="data.frame") {
    bed = ref
    ref = readBED(ref)
    if (filter=="promoter")
      i <- which(isEnh(bed[,4]))
    if (filter=="enhancer")
      i <- which(!isEnh(bed[,4]))
  } else {
    if (filter=="promoter")
      i <- which(isEnh(ref$id))
    if (filter=="enhancer")
      i <- which(!isEnh(ref$id))
  }
  if (is.null(group))
    region <- combineRgn(files)
  else
    region <- combineGrpRgn(files, group)
  unique(c(i, which(ref$id %in% filterRgn(region, ref)$id)))
}

combineRgn <- function(files) {
  region <- readBED(files[1])
  if (length(files)>1){
    for(i in 2:length(files))
      region <- union(region, readBED(files[i]))
  }
  region
}

combineGrpRgn <- function(files, grp) {
  intersect(combineRgn(files[grp==unique(grp)[1]]), combineRgn(files[grp==unique(grp)[2]]))
}

filterRgn <- function(region, ref) {
  s <- findOverlaps(region,ref,select="first")
  ref[unique(s[!is.na(s)])]
}

isEnh <- function(id) {
  (grepl("^GH\\d+",id) & nchar(id)>10) | grepl("^chr\\d+_\\d+_\\d+",id)
}

edgeRank <- function(pg, g, cutoff=0.9, min=1e-7, reverse=TRUE) {
  if (is.character(g))
    g <- graph_from_data_frame(read.table(g, stringsAsFactors = FALSE), directed = FALSE)
  else if (class(g)=="data.frame")
    g <- graph_from_data_frame(g, directed = FALSE)
  m <- as.matrix(get.edgelist(g))
  if(reverse) pg=rev(pg)
  m <- matrix(match(m,pg),ncol=2, byrow = FALSE)
  w <- rowSums(m)
  w[is.na(w)]=min
  w[w<quantile(w, cutoff)]=min
  g <- set_edge_attr(g, "weight", value = w)
  setGRank(pg, g, "weight", min)
}

setGRank <- function(pg, g, name, min=1e-7) {
  v <- match(V(g)$name, pg)
  v[is.na(v)]=min
  set_vertex_attr(g, name, value = v)
}

maxG <- function(g1, cutoff=30) {
  unlist(lapply(g1, function(g) vcount(g) > cutoff))
}

minG <- function(g1, cutoff=30) {
  unlist(lapply(g1, function(g) vcount(g) <=cutoff))
}

exportMultinets <- function(g, n=5, steps=4, rewire=FALSE, simplify=TRUE) {
  if (rewire) {
    # rewire vertices with constant probability
    g <- rewire(g, with=each_edge(0.5))
    # shuffle initial weights and assign them randomly to edges
    E(g)$weight <- sample(E(g)$weight)
  }
  if (simplify) g <- igraph::simplify(g)
  gs <- .walktrapClustering(g, n, steps)
  gs2<- lapply(gs[1:3], function(g2) .walktrapClustering(g2,0,steps))
  gs <- c(gs2[[1]][maxG(gs2[[1]])],gs2[[2]][maxG(gs2[[2]])],gs2[[3]][maxG(gs2[[3]])],gs[-c(1:3)])
  gs[order(unlist(lapply(gs,function(v) mean(V(v)$weight))), decreasing = TRUE)]
}

.walktrapClustering <- function(g, n, steps) {
  g1 <- walktrap.community(g, steps=steps)
  tb <- as.data.frame(table(g1$membership))
  id <- if (n>0)
    as.numeric(tb[order(tb$Freq,decreasing=TRUE),][seq_len(n),1])
  else
    as.numeric(tb[order(tb$Freq,decreasing=TRUE),][,1])
  lapply(id, function(i) {
    induced.subgraph(g, vids=which(g1$membership==i))
  })
}

annotNets <- function(gs, sel=c(94,95), n=5) {
  dbs <- listEnrichrDbs()
  annot <- lapply(gs, function(g) enrichr(V(g)$name, dbs[sel,1]))
  lapply(annot, function(x) {
    res <- do.call("rbind",x)
    res <- res[which(is.na(str_locate(res$Term,"PodNet")[,1])),]
    res <- res[which(is.na(str_locate(res$Term,"PluriNet")[,1])),]
    res[order(res$Combined.Score, decreasing = TRUE),]
  })
}

plotMultinets <- function(gs, cols=NULL, vertex.size=10, plot.label=FALSE, interactive=TRUE) {
  if (is.null(cols))
    cols = c("#FFFFFF", "#FFEEEE", "#FFDDDD", "#FFCCCC", "#FFBBBB", "#FFAAAA", "#FF9999", "#FF8888", "#FF7777", "#FF6666", "#FF5555", "#FF4444", "#FF3333", "#FF2222", "#FF1111", "#FF0000")
  invisible(lapply(gs, function(g) {
    w <- V(g)$weight
    w <- (w - min(w)) / (max(w) - min(w))
    labels <- if(plot.label) {
      V(g)$name
    } else {
      NA
    }
    plot(g, vertex.color=cols[round((length(cols)-1)*w)+1], vertex.size=vertex.size, layout=layout.fruchterman.reingold, vertex.label=labels)
    .interactive(interactive)
  }))
}

exportPCs <- function (res, file) {
  gr=split(res$gr,res$gr$seqnames)
  gr=do.call("rbind",lapply(gr,function(x) x[order(x$start),]))
  write.csv(gr,file=file,row.names=FALSE,quote=FALSE)
}

vJoin <- function (v) {
  v <- ifelse(grepl("\n$",v),substr(v,0,nchar(v)-1),v)
  ifelse(grepl(',$',v),substr(v,0,nchar(v)-1),v)
}

exportJSONnets <- function(gs, prefix) {
  s <- '{"subnetworks":['
  for(g in gs) {
    e <- '"edges": ['
    v <- '{"nodes": ['
    edges <- as_edgelist(g)
    if (nrow(edges)>0) {
      for(j in 1:nrow(edges)){
        e <- paste0(e,'{"source":"',edges[j,1],'","target":"',edges[j,2],'","networks":["Network"]},',"\n")
      }
    }
    vn <- V(g)$name
    for(i in seq_len(length(vn))) {
      v <- paste0(v,'{"name":"',vn[i],"\"},\n")
    }
    s <- paste0(s,vJoin(v),'],',vJoin(e),']},')
  }
  write(paste0(vJoin(s),']}'), file=paste0(prefix,"nets.json"))
}

getJSONarray <- function (df, field, n=5) {
  paste0(df[1:n,field],collapse='","')
}

getJSONpathways <- function (df, n=5) {
  x <- str_split(df[1:n,'Genes'],';')
  df <- data.frame(gene=unlist(x),id=unlist(lapply(1:length(x),function(i) rep(i-1,length(x[[i]])))))
  x <- lapply(split(df, df$gene), function(s) paste0(s$id,collapse=','))
  df <- data.frame(gene=unlist(names(x)),id=unlist(x))
  paste0(apply(df,1,function(s) paste0('"',s[1],'":[',s[2],']')),collapse=',')
}

getJSONgenes <- function (df, n=5) {
  paste0(str_replace_all(df[1:n, "Genes"], ';', '","'),collapse='"],["')
}

exportJSONpathways <- function(enriched, prefix, n=5) {
  s <- '['
  for(df in enriched) {
    n <- min(n, nrow(df))
    s <- paste0(s,'{"pathways": ["',getJSONarray(df,'Term',n),'"],"pvals": ["',getJSONarray(df,'Adjusted.P.value',n),'"],"genes": [["',getJSONgenes(df,n),'"]]},')
  }
  write(paste0(vJoin(s),']'), file=paste0(prefix,"pathways.json"))
}

getId <- function(id, type = c("gene", "all", "enhancer", "promoter", "original"), replacement=''){
  type = match.arg(type)
  if (type == "gene") {
    id=id[!isEnh(id)]
    id=unique(gsub("_\\d+",replacement,id))
  } else if (type == "all") {
    id=unique(gsub("_\\d+",replacement,id))
  } else if (type == "enhancer") {
    id=id[isEnh(id)]
  } else if (type == "promoter") {
    id=id[!isEnh(id)]
  }
  id[id!='0']
}

getPromId <- function(gr, pc="PC1"){
  getId(gr[order(abs(gr[,pc]),decreasing=TRUE),'id'], type="gene")
}

.means <- function(x){
  if (length(dim(x)) < 2L)
    x
  else
    rowMeans(x)
}

plotD <- function(Dobs, labels, captions=NULL, title="", saveFile=FALSE){
  if (saveFile) png(paste0(title,'.png'))
  captions <- if (is.null(captions))
    rep("",ncol(Dobs))
  else
    paste0(" (",captions,")")
  par(mfrow=c(2,4))
  for(i in 1:ncol(Dobs))
  barplot(prcomp(Dobs)$rotation[,i], names=labels, main=paste0("PC",i,captions[i]), las=2)
  par(mfrow=c(1,1))
  if (saveFile) dev.off()
}

exportD <- function(Dobs, labels, prefix){
	x <- t(prcomp(Dobs)$rotation)
	colnames(x) <- labels
	write.csv(x,file=paste0(prefix,"comp.csv"),row.names=F,quote=F)
}

.trapz <- function (x, y) {
    x[length(x)] = 1
    idx = 2:length(x)
    as.double((x[idx] - x[idx - 1]) %*% (y[idx] + y[idx - 1]))/2
}

posId <- function(x, pos, sort=TRUE){
  x <- match(pos,x)
  names(x) <- pos
  if (sort)
  sort(x[!is.na(x)], decreasing = FALSE)
  else
  x[!is.na(x)]
}

#this is an accessory function, that determines the number of points below a diagnoal passing through [x,yPt]
numPts_below_line <- function(myVector,slope,x){
	yPt <- myVector[x]
	b <- yPt-(slope*x)
	xPts <- 1:length(myVector)
	return(sum(myVector<=(xPts*slope+b)))
}

calcCutoff <- function(inputVector){
 	inputVector <- sort(inputVector)
	inputVector[inputVector<0]<-0 #set those regions with more control than ranking equal to zero
	slope <- (max(inputVector)-min(inputVector))/length(inputVector) #This is the slope of the line we want to slide. This is the diagonal.
	xPt <- floor(optimize(numPts_below_line,lower=1,upper=length(inputVector),myVector= inputVector,slope=slope)$minimum) #Find the x-axis point where a line passing through that point has the minimum number of points below it. (ie. tangent)
	y_cutoff <- inputVector[xPt] #The y-value at this x point. This is our cutoff.
}

getAUC <- function (x, pos){
    len <- length(x)
    fx <- posId(x, pos)
    ex <- .ecdf(fx, len)
    s <- seq(0, len, 1)
    1 - .trapz(ex(s), s/len)
}

plotRank <- function(x, pos, title="", file=NA, perc=1, captions=NULL){
  cols <- colorRampPalette(c('red','blue','green'))(length(x))
  len <- length(x[[1]])
  fx <- posId(x[[1]],pos)
  auc <- rep(0, length(x))
  if (!is.na(file)) png(file)
  ex <- .ecdf(fx, len)
  plot(ex, verticals = TRUE, do.points = FALSE, col=cols[1], main=title, cex=2, lwd=2)
  L <- round(perc*len)
  s <- seq(0,L,1)
  auc[1] <- 1 - .trapz(ex(s), s/L)
  if (length(x)>1){
    for(i in 2:length(x)) {
      fx <- posId(x[[i]][x[[i]] %in% x[[1]]],pos)
      ex <- .ecdf(fx, len)
      lines(ex, verticals = TRUE, do.points = FALSE, col=cols[i], lwd=2)
      auc[i] <- 1 - .trapz(ex(s), s/L)
    }
  }
  if (is.null(captions)) captions <- names(x)
  legend("bottomright", legend=paste(captions,paste(round(100*auc, 2), "%", sep="")), col=cols, lty=1, cex=1.4)
  if (!is.na(file)) dev.off()
  names(auc) <- captions
  auc
}

writeRank <- function(id, id2, prefix, header=c("PromEnh","PromOnly")) {
  x <- cbind(id,id2)
  colnames(x) <- header
  write.csv(x, paste0(prefix,"rank.csv"), quote=FALSE, row.names=FALSE)
}

writeMarkers <- function(m, prefix) {
  x <- do.call("rbind",lapply(seq_len(length(m)),function(i) cbind(m[[i]],names(m)[i])))
  colnames(x) <- c("name","source")
  write.csv(x, paste0(prefix,"marker.csv"), quote=FALSE, row.names=FALSE)
}

.lw <- function (x) length(which(x))

.ecdf <- function (x, n) {
    x <- sort(x)
    if (n < 1)
        stop("'x' must have 1 or more non-missing values")
    vals <- unique(x)
    rval <- stats::approxfun(c(vals,n), c(cumsum(tabulate(match(x, vals)))/length(x),1),
        method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
    class(rval) <- c("ecdf", "stepfun", class(rval))
    assign("nobs", n, envir = environment(rval))
    attr(rval, "call") <- sys.call()
    rval
}

.smallPositiveNum <- function(len, x=1e-99) {
  runif(len, min = x, max = x*2)
}

.minNonZero <- function(x) {
  x[x==0] <- min(x[x>0])
  x
}

pageRank <- function(pcs, x, damping = 0.85, dWeight=1e-99, fun = NULL, maxTry=10, gene.id=TRUE, rewire=FALSE, statLog=NULL){
  id <- as.character(pcs[,'id'])
  promId <- getId(id, type="promoter")
  if (!all(promId %in% id))
    stop("Promoter IDs should be a subset of input IDs.")
  if (class(x)=="character")
    x <- read.table(x, stringsAsFactors = FALSE)
  if (class(x)=="data.frame") {
    x <- x[x[,1] %in% id & x[,2] %in% id,]
    ## averaging over enhancer density
    # x[,3] <- x[,3]/table(x[,2])[x[,2]]
    colnames(x) <- c("a","b")
    promId <- promId[!(promId %in% c(x[,1],x[,2]))]
    ## add a dummy node 0 and unappeared gene id into network and assign weights for edges start with 0 small values
    d <- rep(0,length(promId))
    d1<- .smallPositiveNum(length(promId), dWeight)
    w <- c(x[,3],d1)
    w[which(w==0)] = .smallPositiveNum(.lw(w==0), dWeight)
    x <- rbind(x[,1:2], data.frame(a=factor(d),b=promId))
    if (max(w) > 1)
      w <- w/max(w)
    g <- graph_from_data_frame(x, directed=TRUE)
  }
  if (!is.null(statLog)) {
    dg1 <- degree(g,mode="in")
    dg2 <- degree(g,mode="out")
    df <- cbind(as.matrix(summary(dg1[!isEnh(names(dg1))])), as.matrix(summary(dg2[isEnh(names(dg2))])))
    colnames(df) <- c("prom","enh")
    write.csv(df, file=paste0(statLog,"stat.csv"), quote=FALSE, row.names=FALSE)
  }
  pcd <- grep('PCx|PC\\d+',names(pcs), ignore.case = TRUE)
  vs <- vector(mode = "list", length = length(pcd))
  names(vs) <- names(pcs)[pcd]
  for(pc in pcd) {
  wt <- abs(pcs[,pc])
  id <- id[order(wt, decreasing = FALSE)]
  ix <- id %in% V(g)$name
  id <- id[ix]
  id <- c('0',id)
  vr <- match(V(g)$name, id)
  vr[is.na(vr)] <- 1
  dt <- c(dWeight, sort(wt, decreasing = FALSE)[ix])
  dt[dt==0] <- dWeight
  if (rewire) {
    # rewire vertices with constant probability
    g <- rewire(g, with=each_edge(prob=0.5))
    # shuffle initial weights and assign them randomly to edges
    #E(g)$weight <- sample(E(g)$weight)
  }
  iTry <- 0
  while(TRUE){
    v <- try(page_rank(g, algo="arpack", personalized=dt[vr], weights=w, directed = TRUE, damping = damping), silent = TRUE)
    iTry <- iTry + 1
    if(class(v) != "try-error" | iTry>maxTry) break
    w[which.min(w)] <- .smallPositiveNum(1)
  }
  v <- v$vector
  v <- v[order(v, decreasing = TRUE)]
  message(paste0(.lw(v>calcCutoff(v[!isEnh(names(v))]))," promoters are highly ranked in ",names(pcs)[pc]))
  vs[[names(pcs)[pc]]] <- if (gene.id) {
    getId(names(v), type="gene")
  } else {
    v <- v[!isEnh(names(v))]
    names(v) <- gsub("_\\d+", "", names(v))
    v
  }
  }
  vs
}

.interactive <- function(interactive) {
  if (interactive)
    invisible(readline(prompt="Press [enter] to continue"))
}

.condense <- function(d) {
  j=1
  for(i in sort(unique(d))) {
    d[d==i]=j
    j=j+1
  }
  d
}

extractData <- function(data, meta, datasets, groups=c(1,2), logTransform=TRUE, normalize=TRUE) {
  j <- meta[,'group'] %in% groups & meta[,'dataset'] %in% datasets
  meta <- meta[j,]
  data <- data[,j]
  i <- order(meta[,'group'])
  data <- data[,i]
  colnames(data) <- meta[i,'file']
  if (normalize)
    data <- normalize.quantiles(data)
  if (logTransform)
    data <- log(data+1)
  round(data)
}

writeJSONList <- function(data, file){
  write(paste0('["',
  paste0(c(t(data)),collapse='","'),
  '"]'),file=file)
}

writeJSONArray <- function(data, file){
  write(paste0('["',
  paste0(data,collapse='","'),
  '"]'),file=file)
}

writeData <- function(data, groupLabels, prefix, n=61, width=2000, name=NULL, intTemp=FALSE){
  filename <- prefix
  data <- order2(data)
  code <- unlist(strsplit("0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz",""))
  if (length(code) < n)
    stop("do not support so many colors")
  name <- if(is.null(name))
    colnames(data)
  else
    name
  i <- seq(which(grepl('PCx',name))+1,length(name))
  name2<- paste(str_extract(name,"WGB|fractional|Bisulfite|H3K\\w+"),groupLabels[as.integer(str_extract(name,"^\\d+"))], sep = '-')[i]
  name <- paste(groupLabels[as.integer(str_extract(name,"^\\d+"))],str_extract(name,"WGB|fractional|Bisulfite|H3K\\w+"), sep = '-')[i]
  name2<- str_replace(name2,"WGB|fractional|Bisulfite","Methylation")
  name <- str_replace(name, "WGB|fractional|Bisulfite","Methylation")
  ix   <- order(name2)
  gr   <- data[,1:3]
  id   <- data[,'id']
  dpc  <- data[,grepl("PC",colnames(data))]
  data <- data[,i][,ix]
  k <- n/(max(data)-min(data))
  data <- apply(round((data-min(data))*k), 2, function(i) code[i+1])
  data <- apply(data, 1, function(i) paste0(i,collapse=""))
  write.csv(name[ix], file=paste0(filename,'name.csv'), quote=FALSE, row.names = FALSE)
  write.csv(cbind(id=overlapWins(gr, width), gr, name=id, dpc, data), file=paste0(filename,'seq.csv'), quote=FALSE, row.names = FALSE)
  if (intTemp){
  writeLines('{"Prom":[],"Enh":[],"DefaultType":0,"Type":["no capture Hi-C data available"]}',con=paste0(filename,'int.json'))
  writeLines('P,E,T',con=paste0(filename,'int.csv'))}
}

overlapWins <- function(gr, width=2000) {
  findOverlaps(makeGRangesFromDataFrame(gr),makeBedWins(gr, width=width),select="first")
}

makeBedWins <- function(bed, width=2000) {
  bed[,1] <- as.character(bed[,1])
  bed <- bed[!grepl('_',bed[,1]),]
  df <- data.frame(matrix(unlist(lapply(split(bed,bed[,1]), function(x){
    c(as.character(x[1,1]),1,max(x[,3]))
  })), ncol = 3, byrow = TRUE), stringsAsFactors = FALSE)
  if (any(is.na(df[, 1])))
    df <- df[-which(is.na(df[,1])),]
  colnames(df) <- c("seqnames","start","end")
  df$start <- as.integer(df$start)
  df$end <- as.integer(df$end)
  unlist(tile(makeGRangesFromDataFrame(df), width=width))
}

.any <- function(x, n) .lw(x)>=n

dPCA <- function(meta, bed, data, sampleId=NULL, groups=1:2, datasets=NULL, transform=NULL, normlen=NULL, minlen=50, lambda=0.15, fun=function(x) sqrt(mean(x^2)), datasetLabels=NULL, groupLabels=NULL, qnormalize=TRUE, qnormalizeFirst=FALSE, normalize=FALSE, verbose=FALSE, interactive=FALSE, useSVD=FALSE, saveFile=FALSE, processedData=TRUE, removeLowCoverageChIPseq=FALSE, removeLowCoverageChIPseqProbs=0.1, dPCsigns=NULL, nPaired=0, nTransform=0, nColMeanCent=0, nColStand=0, nMColMeanCent=1, nMColStand=0, dSNRCut=5, nUsedPCAZ=0, nUseRB=0, dPeakFDRCut=0.5) {
  #if (class(bed)=="data.frame")
  #  bed = readBED(bed)
  #data = data[match(bed[,4],rownames(data)),]
  colnames(bed) <- c("seqnames","start","end","id")
  datac <- as.character(meta[,'file'])
  if (!is.null(datasets)) {
    wId <- which(meta[,'dataset'] %in% datasets)
    meta <- meta[wId,]
    data <- data[,wId]
    datac<- datac[wId]
  }
  if (!is.null(groups)) {
    wId <- which(meta[,'group'] %in% groups)
    meta <- meta[wId,]
    data <- data[,wId]
    datac<- datac[wId]
  }
  if (!is.null(sampleId)) {
    meta <- meta[sampleId,]
    data <- data[,sampleId]
    datac<- datac[sampleId]
  }
  len <- abs(bed[,3]-bed[,2])+1
  groupId <- as.integer(meta[,'group'])
  datasetId <- as.integer(meta[,'dataset'])
  datac <- paste(groupId, datac, sep = '-') 
  condId <- unique(groupId)
  nGroupNum <- length(condId)
  nDatasetNum <- length(unique(datasetId))
  nSampleNum <- nrow(meta)
  sampleName <- rep("",nSampleNum)
  repNum <- matrix(unlist(lapply(condId, function(i){
    x <- datasetId[groupId==i]
    unlist(lapply(unique(x), function(j){
      .lw(x==j)
    }))
  })),nrow=nGroupNum, byrow=TRUE)
  if (is.null(datasetLabels)) {
    datasetLabels <- seq_len(nDatasetNum)
  }
  if (is.null(groupLabels)) {
    groupLabels <- seq_len(nGroupNum)
  }
  if (!is.null(minlen)) {
    wId <- which(len>minlen)
    bed <- bed[wId,]
    len <- len[wId]
    data<- data[wId,]
  }
  if (!is.null(normlen)) {
    wId <- which(meta[,'dataset'] %in% normlen)
    tdata <- data*1000 / len
    data[,wId] <- tdata[,wId]
  }
  if (normalize) {
    sft <- colSums(data)
    data <- data %*% diag(max(sft)/sft)
  }
  if (removeLowCoverageChIPseq) {
    qv <- max(1,quantile(data, probs=removeLowCoverageChIPseqProbs))
    qi <- apply(data,1,function(x) .any(x>qv, 2))
  }
  if (qnormalizeFirst) {
    data <- normalize.quantiles(data)
  }
  if (!is.null(transform)) {
    wId <- which(meta[,'dataset'] %in% transform)
    tdata <- apply(data,2,function(tmp){
      #tmp <- tmp + 1e-9
      #bc <- boxcox(do.call("lm",list(tmp~1)), plotit = verbose)
      #lam <- bc$x[which.max(bc$y)]
      #if (verbose) print(lam)
      .boxcoxTrans(tmp, lam=lambda)
    })
    data[,wId] <- tdata[,wId]
  }
  if (removeLowCoverageChIPseq) {
    med <- min(data)
    data[!qi,] <- runif(ncol(data), med-0.00001, med)
  }
  if (!qnormalizeFirst & qnormalize) {
    data <- normalize.quantiles(data)
  }
  if (verbose) {
    boxplot(data)
  }
  #bed <- as.data.frame(bed)
  nLociNum <- nrow(bed)
  if (useSVD) {
	g  <- sort(unique(meta$group))
	mx <- sort(unique(meta$dataset))
	d1 <- data[,meta$group==g[1]]
	d2 <- data[,meta$group==g[2]]
	m1 <- meta[meta$group==g[1],'dataset']
	m2 <- meta[meta$group==g[2],'dataset']
	d1 <- do.call("cbind",lapply(mx, function(i) rowMeans(d1[,m1==i])))
	d2 <- do.call("cbind",lapply(mx, function(i) rowMeans(d2[,m2==i])))
	d  <- d1-d2
	d  <- list('PC'=svd(d)$u, 'Dobs'=d, 'proj'=t(data.frame('eigenvalue'=rep(0,ncol(d)),'percent_var'=rep(0,ncol(d)),'SNR'=rep(0,ncol(d)),'sig_prc'=rep(0,ncol(d)))))
  }else{
  d <- dPCA_main_impl(nGroupNum, nDatasetNum, nSampleNum, nPaired, nLociNum, .condense(groupId), .condense(datasetId), repNum, sampleName, bed[,1], bed[,2], bed[,3], data, nTransform, nColMeanCent, nColStand, nMColMeanCent, nMColStand, dSNRCut, nUsedPCAZ, nUseRB, dPeakFDRCut)
  d$proj <- matrix(unlist(lapply(1:4,function(i) d[[i]])),nrow=4,byrow=TRUE,dimnames=list(names(d)[1:4]))
  d$Dobs <- matrix(d$Dobs, ncol=nDatasetNum, byrow=TRUE)
  }
  pc <- matrix(d$PC, ncol=nDatasetNum, byrow=TRUE)
  pcc <- paste0("PC",seq_len(nDatasetNum))
  colnames(pc) <- pcc
  bed <- cbind(bed, pc)
  if (!is.null(dPCsigns))
    pc <- t(t(pc)*dPCsigns)
  bed$PCx <- apply(pc,1,fun)
  colnames(data) <- datac
  if (processedData) bed <- cbind(bed,data)
  list('gr'=bed[order(abs(bed$PC1),decreasing=TRUE),], 'Dobs'=d$Dobs, 'proj'=d$proj)
}

order2 <- function(bed) do.call("rbind",lapply(split(bed,bed[,1]), function(x) x[order(x[,2]),]))

.boxcoxTrans <- function(x, lam) {
  if (lam > 0)
    (x^lam - 1)/lam
  else
    log(x)
}

.findOverlapIndex <- function(gr1, gr2) {
  df <- as.data.frame(findOverlaps(gr1, gr2))
  df <- split(df, df[,1])
  lapply(df, function(x) x[,2])
}

.tblNode <- function(app, name, label="view") {
	paste0('<td><a href="',app,'.html?sample=',name,'" target="_blank">',label,'</a></td>')
}

writeIndexHtml <- function(name, exdir = ".") {
writeLines(paste0('<!DOCTYPE html><html><head><title>',name,'</title><link rel="stylesheet" type="text/css" href="styles/main.css"></head><body><h1>',name,'</h1><table id="rbl"><tr><th>Sample name</th><th>Network enrichment</th><th>Epigenome heatmap</th><th>Rank lists</th><th>Hallmark ROC</th></tr>',
"<tr>",.tblNode("stat",name,name),paste0(sapply(c("net","browse","rank","roc"), function(d) .tblNode(d, name)),collapse=''),"</tr></table></body></html>",collapse=''), con=paste0(exdir,'/index.html'))
}

exportApps <- function(name, markers=NULL, exdir = ".") {
	untar(system.file("data", "html.tar.gz", package="irene"), exdir = exdir)
	writeIndexHtml(name, exdir)
	if (!is.null(markers)) {
		writeMarkers(markers, name)
	}

}

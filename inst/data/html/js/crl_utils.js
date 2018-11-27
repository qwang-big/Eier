function makelist(data){
    var j = 0
		for(var i in data) {
      if (srclist.indexOf(data[i].source) == -1) {
        srclist.push(data[i].source)
        genesrc[j] = []
        j ++
      }
      genesrc[j-1].push(data[i].name)
    }
    for(var i in srclist) {
      d3.select("#lists")
      .append("a")
      .attr("href", "javascript:void(0)")
      .attr("onclick", "highlight("+i+")")
      .text(srclist[i]).append("br")
    }
}

function rank2p(d){
	return shortd(100-Math.max(1,d)/179)
}

function highlight(i){
	chart.selectAll(".bar").remove()
	highlights = genesrc[i]
  plotPt(rank)
}

function highlightgene(data){
  var bar = chart
    .append("g")
    .datum(data)
    .attr("class", "bar")
    .append("svg:rect")
    .on('click', function(d,i) {showPtinfo(d.name)})
    .attr("x", function(d) { return x(d.rank) })
    .attr("y", function(d) {return y(d.rank2)})      
    .attr("height", 9)
    .attr("width", 9)
    .attr("fill", "#f00")
  .on("mouseover", function(d){
	tooltip.html(d.name)
	return tooltip.style("visibility", "visible");
      })
      .on("mousemove", function(){
           return tooltip.style("top", (d3.event.pageY-10)+"px").style("left",(d3.event.pageX+10)+"px")
      })
      .on("mouseout", function(){
           return tooltip.style("visibility", "hidden")
      })
}

function showPtinfo(d){
	var root = d3.select("#info0")
	root.selectAll("*").remove();
	root.append("span").style("white-space", "nowrap").text(d)
	var table = root.append("table")
	table.append("thead").append("tr")
	    .selectAll("th")
	    .data(["Entity","Distance(bp)"].concat(vpclabels()))
	    .enter().append("th")
	    .append("a").attr("href",function(d,j) { return (j<2)?'javascript:void(0)':'javascript:plotD(0,'+(j-2)+')'})
	    .text(function(d) {
		return d;
	    });
	table.append("tbody").selectAll('tr')
	  .data(d in dPCs ? dPCs[d] : [])
	  .enter()
	  .append("tr")
	  .html(function(d,i) { return tdfill(d, i) })
	legend2(col2, root.append("span")
	  .style("white-space", "nowrap"))
}

function plotPt(data){
  var bar = chart.selectAll(".bar")
    .data(data)
    .enter().append("g")
    .attr("class", "bar");

  bar.append("svg:rect")
    .on('click', function(d,i) {showPtinfo(d.name)})
    .attr("x", function(d) { return x(d.rank) })
    .attr("y", function(d) {return y(d.rank2)})      
    .attr("height", function(d) { return highlights.indexOf(d.name)<0 ? 3 : 9})
    .attr("width", function(d) { return highlights.indexOf(d.name)<0 ? 3 : 9})
    .attr("fill",function(d) { return highlights.indexOf(d.name)<0 ? "#ccc" : "#f00"})
  .on("mouseover", function(d){
	tooltip.html(d.name)
	return tooltip.style("visibility", "visible");
      })
      .on("mousemove", function(){
           return tooltip.style("top", (d3.event.pageY-10)+"px").style("left",(d3.event.pageX+10)+"px")
      })
      .on("mouseout", function(){
           return tooltip.style("visibility", "hidden")
      })
}

function rep(c, n) {
	ret = ''
	for(i=0; i<n; i++)
		ret += c
	return ret
}
function explode(str) {
	ret = ''
	for(i=0; i<str.length; i++) {
		c = str.charAt(i)
		ret += (c>'9') ? c : rep(c1, c-1)
		c1 = c
	}
	return ret
}
function strsplit(x, s, part) {
    if (x.indexOf(s)<0)
      return x
    else if (part==1)
      return x.substr(0,x.indexOf(s))
    else if (part==2)
      return x.substr(x.indexOf(s)+1,x.length)
    else
      return x
  }

function firstd(d) {
	for(var i=0; i<d.length; i++){
	  if (d.charAt(i)=='0'||d.charAt(i)=='.') continue;
	  return (i>d.indexOf(".")?i:d.indexOf("."))+2
	}
}

function shortp(d) {
	return d.indexOf(".")>0 ? d.slice(0,firstd(d)) : d
}

function shorte(d) {
	return d.indexOf("e")>0 ? d.slice(d.indexOf("e")) : ""
}

function shortd(d) {
	return shortp(""+d) + shorte(""+d)
}

function ashortd(d) {
	d = shortd(d)
	return d.charAt(0)=="-" ? d.slice(1) : d
}

function getReadableIndex(d, w) {
	return (d/1000000).toFixed(2)+"Mb"
}

function zoomIn(d) {
	scale=Math.min(20, scale*2)
	heatmapImpl()
	return false
}

function zoomOut(d) {
	scale=Math.max( 1, scale/2)
	heatmapImpl()
	return false
}

function heatmap(d) {
  curMap = d
  curMap.s = Math.max(0,d.id-200)
  curMap.e = d.id+200
  d3.select("head").select('title').text(d.name)
  heatmapImpl()
}

function getParameterMap () {
    if (window.location.href.indexOf('?') === (-1)) {
        return {};
    }
    var qparts = window.location.href.split('?')[1].split('&'),
        qmap   = {};
    qparts.map(function (part) {
        var kvPair = part.split('='),
            key    = decodeURIComponent(kvPair[0]),
            value  = kvPair[1];

        //handle params that lack a value: e.g. &delayed=
        qmap[key] = (!value) ? '' : decodeURIComponent(value);
    });
    return qmap;
}

function legend(cols, id, left, right) {
	var panel = d3.select(id)
  	panel.selectAll("*").remove()
	panel.append('span').text(left)
	var svg = panel.append('svg').attr("width", 5*cols.length).attr("height", 20)
	for(i in cols)
	svg.append('rect')
	.attr('x',5*i)
	.attr('y',0)
	.attr('width',5)
	.attr('height',20)
	.attr('fill', cols[i] )
	panel.append('span').text(right)
}

function legend2(cols, panel) {
	panel.append('span').text(shortd(settings.low))
	var svg = panel.append('svg').attr("width", 5*col2.length).attr("height", 20)
	for(i in cols)
	svg.append('rect')
	.attr('x',5*i)
	.attr('y',0)
	.attr('width',5)
	.attr('height',20)
	.attr('fill', col2[i] )
	panel.append('span').text(shortd(settings.high))
}

function round2p(num) {    
    return Math.round(num * 100);
}

function getscales(seqs, s, e) {
	var a=code.length-1, b=0
  for(var i=s; i<e; i++) {
		if (i in seqs){
			x = getscale(seqs[i])
			a = Math.min(a,x[0])
			b = Math.max(b,x[1])
		}
	}
	return [a,Math.floor((a+b)/2),b]
}

function getscale(seq) {
	var a=code.length-1, b=0
	for(var i=0; i<seq.length; i++) {
		x = code.indexOf(seq.charAt(i))
		a = Math.min(a,x)
		b = Math.max(b,x)
  }
	return [a,b]
}

function getIntEnh(genename, ti) {
	var res = []
	if (genename in inta && ti in inta[genename])
	for(var i in inta[genename][ti]){
		res.push(intc.Enh[inta[genename][ti][i]])
	}
	return res
}

function redirecturl(id) {
	return isEnh(id) ? 'redirect.html?ghid='+id : "http://www.genecards.org/Search/Keyword?queryString="+id
}

function heatmapImpl() {
	  var d = curMap
	  var root = d3.select("#heatmap")
	  root.selectAll("*").remove()
    var tooltip = d3.select("#tooltip")
    color = d3.scale.linear().domain([d3.min(dPC),0,d3.max(dPC)]).range(["#00f", "#fff", "#f00"]).nice()
	  legend(cols, "#legend", "Epigenome intensity percentage scale: 0%", "100%")
	  legend(col2, "#legend2", "PC1: "+shortd(d3.min(dPC)), shortd(d3.max(dPC)))
	  var h=20, y=0, x=0, W=1500, H=h*names.length
	  var s = parseInt(d.s)
	  var e = parseInt(d.e)
	  var p = e - s, n = parseInt((p - p/scale)/2)
	  s += n
	  e -= n
	  p  = p/scale
	  var w = Math.max(W/p, 1)
	  var start=d.start-p*1000, end=d.end+p*1000
	  var plot = root.append('svg').attr("width", W+150).attr("height", H+h*2)
	  plot.append('text')
	  .attr('x',0)
	  .attr('y',h-4)
	  .text(getReadableIndex(start))
	  .style("font-size", "16px")
	  plot.append('text')
	  .attr('x',W)
	  .attr('y',h-4)
	  .text(getReadableIndex(end))
	  .style("font-size", "16px")
    plot.append('rect')
    .attr('x',0)
    .attr('y',h)
    .attr('width',W)
    .attr('height',h)
    .style("z-index", -1)
    .style('fill', '#eee')
    //.style("fill-opacity",0)
    //.style('stroke-width', 0.3)
    //.style('stroke', '#000')
	  for(var i=s; i<e; i++) {
			if (i in TAD)
	      plot.append('a')
        .append('rect')
	      .attr('x',x)
	      .attr('y',h)
	      .attr('width',w)
	      .attr('height',h)
	      .attr('fill', TAD[i]%2==0 ? '#fc0' : '#c0f' )
			x += w
		}
		x = 0
    color2= d3.scale.linear().domain(getscales(data, s, e)).range(["#BEBEBE", "#FFFFFF", "#FF0000"]).nice()
	  for(var i=s; i<e; i++) {
      y = h
			if (i in data){
	      plot.append('a')
		    .attr('xlink:href',redirecturl(strsplit(ids[i],'_',1)))
		    .attr('target',"_blank")
        .append('rect')
        .datum(i)
	      .attr('x',x)
	      .attr('y',y)
	      .attr('width',w)
	      .attr('height',h)
	      .attr('fill', i in dPC ? color(dPC[i]) : '#fff' )
			  .style('stroke-width', 0.1)
			  .style('stroke', '#000')
	      .on("mouseover", function(j){
            tooltip.html(ids[j]+": "+(dPC[j]!== undefined?dPC[j].toFixed(2):"undefined"))
	          return tooltip.style("visibility", "visible");
        })
	      .on("mousemove", function(){
	           return tooltip.style("top", (d3.event.pageY-10)+"px").style("left",(d3.event.pageX+10)+"px")
	      })
	      .on("mouseout", function(){
	           return tooltip.style("visibility", "hidden")
	      })
				y += h
	      for(var j=0; j<data[i].length; j++) {
		      plot.append('rect')
      	  .datum(round2p(code.indexOf(data[i].charAt(j))/code.length))
		      .attr('x',x)
		      .attr('y',y)
		      .attr('width',w)
		      .attr('height',h)
		      .attr('fill', color2(code.indexOf(data[i].charAt(j))))
			    .on("mouseover", function(s){
		          tooltip.html(s+"%")
			        return tooltip.style("visibility", "visible");
		      })
			    .on("mousemove", function(){
			         return tooltip.style("top", (d3.event.pageY-10)+"px").style("left",(d3.event.pageX+10)+"px")
			    })
			    .on("mouseout", function(){
			         return tooltip.style("visibility", "hidden")
			    })
					if (i == d.id)
						plot.append('svg:path')
					  .attr('transform', function(d) { return 'translate(' + (x+w/2) + ',' + (h-5) + ')' })
						.attr("d", d3.svg.symbol()
						.type(function(d) { return d3.svg.symbolTypes[4] })
						.size(function(d) { return 40 }))
						.attr('fill', function(d) { return "#f00" })
					else if (intc != undefined && getIntEnh(d.name,intc.sel).indexOf(ids[i])!=-1)
						plot.append('svg:path')
					  .attr('transform', function(d) { return 'translate(' + (x+w/2) + ',' + (h-5) + ')' })
						.attr("d", d3.svg.symbol()
						.type(function(d) { return d3.svg.symbolTypes[4] })
						.size(function(d) { return 40 }))
						.attr('fill', function(d) { return "#A020F0" })
				  y += h
		    }
      }else{
	      plot.append('rect')
	      .attr('x',x)
	      .attr('y',h+h)
	      .attr('width',w)
	      .attr('height',H)
	      .attr('fill', '#C7EDCC')
      }
      x += w
		}
    y = h+h
	  plot.append('text')
	  .attr('x',x+w)
	  .attr('y',y-9)
	  .text("promoter/enhancer")
	  .style("font-size", "16px")
    var group = strsplit(names[0],'-',1)
    var dataset = strsplit(names[0],'-',2)
    for(var j=0; j<names.length; j++) {
		//if (typeof jsonUrl!="undefined" && jsonUrl)
		  plot.append('a')
		  .attr('xlink:href',epggUrl+jsonUrl+"&coordinate="+d.chr+":"+start+"-"+end)
		  .attr('target',"_blank")
		  .append('text').attr('x',x+w)
		  .attr('y',11+y)
		  .text(names[j])
		  .style("font-size", "16px")
      if (group !== strsplit(names[j],'-',1)){
        plot.append('line')
	      .attr('x1',0)
	      .attr('y1',y)
	      .attr('x2',W)
	      .attr('y2',y)
			  .style('stroke-width', 0.1)
			  .style('stroke', '#000')
      }
      group = strsplit(names[j],'-',1)
      if (dataset !== strsplit(names[j],'-',2)){
        plot.append('line')
	      .attr('x1',0)
	      .attr('y1',y)
	      .attr('x2',W)
	      .attr('y2',y)
			  .style('stroke-width', 0.3)
			  .style('stroke', '#000')
      }
      dataset = strsplit(names[j],'-',2)
			y += h
    }
  }

function tdfill(d, i) {
	var colord = d3.scale.linear()
	  .domain([settings.low, 0, settings.high])
	  .range(["#00f", "#fff", "#f00"])
	  .nice();
	var s = isEnh(d.name) ? 'redirect.html?ghid=' : 'browse.html?sample='+sample+'&gene=';
	var str = '<td><a href="'+s+strsplit(d.name,'_',1)+'" target="_blank">'+d.name+'</a></td><td>'+d.dist+'</td>'
	for(var j=1; j<=settings.pc; j++)
		str += '<td style="background:'+colord(parseFloat(d["PC"+j]))+'">&nbsp;</td>'
	return str
}

function trimstr(str, l) {
	return str.substr(0,str.length-l)
}

function vpclabels() {
	var foo = [];
	for (var i = 1; i <= settings.pc; i++) {
	   foo.push('PC'+i);
	}
	return foo
}

function cloned(obj, name, d) {
	var target = {};
	for (var i in obj) {
		if (!Object.prototype.hasOwnProperty.call(obj, i)) continue;
		target[i] = obj[i]
	}
	target[name] = d
	return target
}

function range(d, i, n) {
	var a = []
	for(j=Math.max(d-i,0); j<Math.min(d+i, n); j++)
		a.push(j)
	return a;
}

function showtooltip(text) {
	tooltip.html(text)
}

function plotD(i, j){
	var w = 200, xm = 260, y=16
	var vec = eigenVec[j]
	var keys = Object.keys(vec)
	var svg = d3.select("#info"+i).append('svg').attr("width", xm*2).attr("height", 20*keys.length)
	svg.append('text')
	.attr('x',1)
	.attr('y',12)
	.text("PC"+(j+1))
	.style("font-size", "16px")
	svg.append('line')
    .attr('x1',xm-w)
    .attr('y1',12)
    .attr('x2',xm+w)
    .attr('y2',12)
	  .style('stroke-width', 0.3)
	  .style('stroke', '#000')
	svg.append('line')
    .attr('x1',xm-w)
    .attr('y1',7)
    .attr('x2',xm-w)
    .attr('y2',12)
	  .style('stroke-width', 0.3)
	  .style('stroke', '#000')
	svg.append('line')
    .attr('x1',xm+w)
    .attr('y1',7)
    .attr('x2',xm+w)
    .attr('y2',12)
	  .style('stroke-width', 0.3)
	  .style('stroke', '#000')
	svg.append('line')
    .attr('x1',xm)
    .attr('y1',7)
    .attr('x2',xm)
    .attr('y2',12)
	  .style('stroke-width', 0.3)
	  .style('stroke', '#000')
  svg.append('text')
		.attr('x',xm-w-10)
		.attr('y',11)
		.text("-1")
		.style("font-size", "12px")
  svg.append('text')
		.attr('x',xm+w+1)
		.attr('y',11)
		.text("1")
		.style("font-size", "12px")
  svg.append('text')
		.attr('x',xm+1)
		.attr('y',11)
		.text("0")
		.style("font-size", "12px")
	for(var k in keys) {
		f = parseFloat(vec[keys[k]])
		svg.append('text')
		.attr('x',1)
		.attr('y',10*k+y+12)
		.text(keys[k])
		.style("font-size", "16px")
		svg.append('rect')
		.attr('x',f>0?xm:xm + w*f)
		.attr('y',10*k+y)
		.attr('width',w*Math.abs(f))
		.attr('height',10)
		.attr('fill', f>0?'#f00':'#00f' )
		y += 5
	}
}

function showinfo(d, i){
	var root = d3.select("#info"+i)
	root.selectAll("*").remove();
	root.append("span").style("white-space", "nowrap").text(d.name)
	var table = root.append("table")//.attr("id", "tableID")
	table.append("thead").append("tr")
	    .selectAll("th")
	    .data(["Entity","Distance(bp)"].concat(vpclabels()))
	    .enter().append("th")
	    .append("a").attr("href",function(d,j) { return (j<2)?'javascript:void(0)':'javascript:plotD('+i+','+(j-2)+')'})
	    .text(function(d) {
		return d;
	    });
	table.append("tbody").selectAll('tr')
	  .data(d.name in dPCs ? dPCs[d.name] : [])
	  .enter()
	  .append("tr")
	  .html(function(d,i) { return tdfill(d, i) })
	//$('#tableID').dataTable({searching:false, paging:true, ordering:true, info:false, "order": [[ 1, "desc" ]]});
	legend2(col2, root.append("span")
	  .style("height", "25px")
	  .style("white-space", "nowrap")
	  .attr("id", "info-legend"+i))
}

function strokecolor(d) {
    return ("expr" in d) ? dColor(d.expr) : "#000";
}

function getNetstr(d, i) {
	str = ''
	for(var j in d.nodes) {
		if('name' in d.nodes[j] && d.nodes[j].name != undefined)
		str += d.nodes[j].name+"\r\n"
	}
	return str
}

function isEnh(x) {
	return x.indexOf("_")==-1 && x.indexOf("GH")==0
}

function setRanks(data) {
	prom =[], rank =[]
	for(i in data) {
		prom[data[i].PromOnly]=data.length-i
	}
	for(i in data) {
		rank[data[i].PromEnh]={"name": data[i].PromEnh, "rank": data.length-i, "rank2": prom[data[i].PromEnh]}
	}

}

function setRank2(data) {
	prom =[], rank =[]
	for(i in data) {
		prom[data[i].PromOnly]=data.length-i
	}
	for(i in data) {
		rank.push({"name": data[i].PromEnh, "rank": data.length-i, "rank2": prom[data[i].PromEnh]})
	}

}

function setPCs(d, l) {
	var s=[], n=10
	for(i=1;i<100;i++) {
		if(Object.keys(d[0]).includes("PC"+i))
			settings.pc = i
		else
			break
	}
	for (var i in d) {
		d[i].start = parseInt(d[i].start)
		d[i].end = parseInt(d[i].end)
		d[i].PC1 = parseFloat(d[i].PC1)
		if (d[i].name.indexOf('_')>=0) {
			s.push({i: parseInt(i), seqnames: d[i].seqnames, start: d[i].start-l, end: d[i].end+l, id: trimstr(d[i].name,2)})
			if (settings.high < d[i].PC1) settings.high = d[i].PC1
			if (settings.low  > d[i].PC1) settings.low  = d[i].PC1
		}
	}
	for (var i in s) {
		dPCs[s[i].id] = []
		for(var j=Math.max(s[i].i-n,0); j<Math.min(s[i].i+n, d.length-1); j++) {
			if ((d[j].name.indexOf(s[i].id)==0||d[j].name.indexOf('_')<0) && s[i].seqnames==d[j].seqnames && s[i].start<d[j].start && s[i].end>d[j].end) {
				dPCs[s[i].id].push(cloned(d[j],"dist",d[j].start-s[i].start-l))
			}
		}
	}
}


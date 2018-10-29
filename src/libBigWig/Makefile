CC ?= gcc
AR ?= ar
RANLIB ?= ranlib
CFLAGS ?= -g -Wall -O3 -Wsign-compare
LIBS = -lcurl -lm -lz
EXTRA_CFLAGS_PIC = -fPIC
LDFLAGS =
LDLIBS =
INCLUDES = 

all: libBigWig.a

OBJS = io.o bwValues.o bwRead.o bwStats.o bwWrite.o

.c.o:
	$(CC) -I. $(CFLAGS) $(INCLUDES) $(EXTRA_CFLAGS_PIC) -c -o $@ $<

libBigWig.a: $(OBJS)
	-@rm -f $@
	$(AR) -rcs $@ $(OBJS)
	$(RANLIB) $@

clean:
	rm -f *.o libBigWig.a


CC=			gcc
CFLAGS=		-g -Wall -O2 #-fno-inline-functions -fno-inline-functions-called-once
CPPFLAGS=
INCLUDES=	
OBJS=		kthread.o rld0.o sys.o diff.o sub.o unpack.o correct.o dfs.o fet.o \
			ksw.o seq.o mag.o unitig.o bubble.o
PROG=		fermi2
LIBS=		-lm -lz -lpthread

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(CPPFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

fermi2:$(OBJS) main.o
		$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

rld0.o:rld0.h
correct.o:kseq.h
seq.o:kseq.h

clean:
		rm -fr gmon.out *.o ext/*.o a.out $(PROG) *~ *.a *.dSYM session*

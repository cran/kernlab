#
# TRON directory.
#

AR     = ar rcv
RANLIB = ranlib 

FILES = dbreakpt.o dcauchy.o dgpnrm.o dgpstep.o dprsrch.o misc.o \
        dspcg.o dtron.o dtrpcg.o dtrqsol.o dprecond.o

CFLAGS = $(OPTFLAGS) 
FFLAGS = $(OPTFLAGS)

src: $(FILES) 
	$(AR) tron.a $(FILES)  
	$(RANLIB) tron.a

clean:
	- rm -f *.o
	- rm -f *.a
	- rm -f *~

.c.o:
	$(CC) $(CFLAGS) -c $*.c



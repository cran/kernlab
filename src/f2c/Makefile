AR     = ar rcv
RANLIB = ranlib 

HEADERS = f2c.h blas.h blasp.h
FILES = dnrm2.o dasum.o daxpy.o dcopy.o ddot.o dscal.o dsymv.o dgemv.o dtrsv.o\
xerbla.o lsame.o dpotf2.o

CFLAGS = $(OPTFLAGS) 
FFLAGS = $(OPTFLAGS)

blas: $(FILES) $(HEADERS)
	$(AR) f2c.a $(FILES)  
	$(RANLIB) f2c.a

clean:
	- rm -f *.o
	- rm -f *.a
	- rm -f *~

.c.o:
	$(CC) $(CFLAGS) -c $*.c



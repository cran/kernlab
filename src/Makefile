CC = gcc
CXXC = g++
CFLAGS = -Wall -g -O3 -fPIC -I$(R_HOME)/include 

all: tron bsvm.o solvebqp.o svm.o kernlab.so 

kernlab.so:  
	$(CXXC) $(CFLAGS) -shared bsvm.o svm.o solvebqp.o \
dtron/tron.a f2c/f2c.a -o kernlab.so -lm $(LIBR)

bsvm.o: bsvm.cpp bsvm.h
	$(CXXC) $(CFLAGS) -c bsvm.cpp
solvebqp.o: solvebqp.c bsvm.h
	$(CC) $(CFLAGS) -c solvebqp.c
svm.o: svm.cpp svm.h 
	$(CXXC) $(CFLAGS) -c svm.cpp
tron: svm.h
	cd dtron; make OPTFLAGS='$(CFLAGS)' CC='$(CC)';
	cd f2c; make OPTFLAGS='$(CFLAGS)' CC='$(CC)';
clean:
	cd f2c;		make clean
	cd dtron;	make clean
	rm -f *~ bsvm.o solvebqp.o bsvm-train bsvm-predict svm-scale

# Makefile to build fea program
# --- macros
CC=gcc
#CFLAGS=  -O3 -g -I./lib/ 
CFLAGS= -g -I./lib/ -I./ -I/group/hpc-share/gsl/1.15/include
OBJECTS= matrixOperations.o prTimesLhood.o main.o
#LIBS = -L/usr/include/gsl -lgsl -L./lib/ -llapack_LINUX -lblas_LINUX -llibf2c -lm
LIBS = -L/group/hpc-share/gsl/1.15/lib -lgsl -L/group/hpc-share/gsl/1.15/lib -l gslcblas -L./lib/ -llapack_LINUX -lblas_LINUX -llibf2c -lm

# --- targetscc
all: fea
fea:$(OBJECTS) 
	$(CC)  -o fea $(OBJECTS) $(LIBS) 

main.o: ./lib/f2c.h matrixOperations.h prTimesLhood.h
	$(CC) $(CFLAGS) -c main.c

matrixOperations.o: ./lib/f2c.h ./lib/blaswrap.h ./lib/clapack.h 
	$(CC) $(CFLAGS) -L/usr/include/gsl -c matrixOperations.c
       
prTimesLhood.o: matrixOperations.h ./lib/f2c.h  
	$(CC) $(CFLAGS) -L/usr/include/gsl -c prTimesLhood.c


# --- remove binary and executable files
clean:
	rm -f fea $(OBJECTS)

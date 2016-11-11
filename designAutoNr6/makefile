# Makefile to build belham_valley_data_reader_auto program
# --- macros
CC=gcc
CFLAGS= -g -I./lib/ -I./
LIBS= -L./lib/ -llapack_LINUX -lblas_LINUX -llibf2c -lm
# --- targetscc
OBJECTS = matrixOperations.o readInZgrid.o designAutoNr.o

all: designAutoNr
designAutoNr:$(OBJECTS) 
	$(CC) -o designAutoNr $(OBJECTS) $(LIBS)

designAutoNr.o: ./lib/f2c.h matrixOperations.h readInZgrid.h
	$(CC) $(CFLAGS) -c designAutoNr.c

matrixOperations.o: ./lib/f2c.h ./lib/blaswrap.h ./lib/clapack.h 
	$(CC) $(CFLAGS)  -c matrixOperations.c

readInZgrid.o: matrixOperations.h
	$(CC) $(CFLAGS) -c readInZgrid.c


# --- remove binary and executable files
clean:
	rm -f designAutoNr $(OBJECTS)

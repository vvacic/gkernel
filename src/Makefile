CC=g++
DEBUG=-g
FLAGS=-Wall -O2
LIBS=-lm


all: parse_pdb make_kernel


aminoacid.o: aminoacid.h aminoacid.cc
	$(CC) $(DEBUG) $(FLAGS) -c aminoacid.cc

string.o: string.h string.cc
	$(CC) $(DEBUG) $(FLAGS) -c string.cc

pdb.o: pdb.h pdb.cc
	$(CC) $(DEBUG) $(FLAGS) -c pdb.cc

protgraph.o: protgraph.h protgraph.cc
	$(CC) $(DEBUG) $(FLAGS) -c protgraph.cc 

gkernel.o: gkernel.h gkernel.cc
	$(CC) $(DEBUG) $(FLAGS) -c gkernel.cc


parse_pdb: parse_pdb.cc aminoacid.o string.o pdb.o protgraph.o
	$(CC) $(DEBUG) $(FLAGS) parse_pdb.cc protgraph.o pdb.o aminoacid.o string.o -o $@

make_kernel: make_kernel.cc aminoacid.o string.o pdb.o protgraph.o gkernel.o
	$(CC) $(DEBUG) $(FLAGS) make_kernel.cc aminoacid.o string.o pdb.o protgraph.o gkernel.o -o $@


clean:
	rm -rf parse_pdb make_kernel *.o *.*~ *~ core *.dSYM


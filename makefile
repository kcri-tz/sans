CC = g++

#CFLAGS = -DMAX_KMER_SIZE=64 -march=native -lbifrost -pthread -lz -lpthread -std=c++11 -Wall -O3 -mcmodel=medium
CFLAGS = -march=native -lbifrost -pthread -lz -lpthread -std=c++11 -Wall -O3 -mcmodel=medium


all:
	$(CC) src/SANS.cpp $(CFLAGS) -o SANS

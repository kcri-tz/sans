CC = g++ -D K=32 -D N=64 # MAX. K-MER LENGTH, NUMBER OF FILES, PLEASE EDIT HERE

SANS: main.o
	$(CC) -o SANS main.o graph.o kmer32.o kmerXX.o color64.o colorXX.o util.o
	rm -rf obj/; mkdir obj/; mv *.o obj/

main.o: src/main.cpp src/main.h graph.o util.o
	$(CC) -c src/main.cpp

graph.o: src/graph.cpp src/graph.h kmer32.o kmerXX.o color64.o colorXX.o
	$(CC) -c src/graph.cpp

kmer32.o: src/kmer32.cpp src/kmer32.h
	$(CC) -c src/kmer32.cpp

kmerXX.o: src/kmerXX.cpp src/kmerXX.h
	$(CC) -c src/kmerXX.cpp

color64.o: src/color64.cpp src/color64.h
	$(CC) -c src/color64.cpp

colorXX.o: src/colorXX.cpp src/colorXX.h
	$(CC) -c src/colorXX.cpp

util.o: src/util.cpp src/util.h
	$(CC) -c src/util.cpp

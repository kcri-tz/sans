# MAX. K-MER LENGTH, NUMBER OF FILES
CC = g++ -O3 -march=native -DmaxK=32 -DmaxN=64 -std=c++14
BF = -lpthread
## IF DEBUG
# CC = g++ -g -march=native -DmaxK=33 -DmaxN=64 -std=c++14

## IF BIFROST LIBRARY SHOULD BE USED
# CC = g++ -O3 -march=native -DmaxK=64 -DmaxN=64 -DuseBF -std=c++14
# BF = -lbifrost -lpthread

# GZ STREAM LIB
CFLAGS = gcc -O3 -march=native

# Wrap Windows / Unix commands
ifeq ($(OS), Windows_NT)
	TD = obj
	MK = rmdir /s /q obj && mkdir obj
	RM = rmdir /s /q obj
	MV = cmd /C move *.o obj
	CP = cp makefile obj

else
	TD = obj/
	MK = mkdir -p obj/
	RM = rm -rf obj/
	MV = ""mv *.o obj/
	CP = cp makefile obj/makefile
	

endif

ifeq ("$(wildcard $(TD))", "")
    RM = @echo ""
endif

SANS: obj/main.o
	$(CC) -o SANS obj/main.o obj/graph.o obj/spinlockMutex.o obj/kmer.o obj/kmerAminoXX.o obj/kmerAmino12.o obj/color.o obj/util.o obj/translator.o obj/cleanliness.o obj/gzstream.o -lz $(BF)
	$(RM)
	$(MK)

obj/main.o: makefile src/main.cpp src/main.h obj/translator.o obj/graph.o obj/util.o obj/cleanliness.o obj/gzstream.o
	$(CC) -c src/main.cpp
	@$(MV)

obj/graph.o: makefile src/graph.cpp src/graph.h obj/kmer.o obj/kmerAmino12.o obj/kmerAminoXX.o obj/color.o obj/spinlockMutex.o
	$(CC) -c src/graph.cpp
	@$(MV)

obj/spinlockMutex.o: src/spinlockMutex.cpp src/spinlockMutex.h
	$(CC) -c src/spinlockMutex.cpp
	@$(MV)

obj/kmer.o: src/kmer.cpp src/kmer.h
	$(CC) -c src/kmer.cpp
	@$(MV)

obj/kmerAmino12.o: src/kmerAmino12.cpp src/kmerAmino12.h obj/util.o
	$(CC) -c src/kmerAmino12.cpp
	@$(MV)

obj/kmerAminoXX.o: makefile src/kmerAminoXX.cpp src/kmerAminoXX.h obj/util.o
	$(CC) -c src/kmerAminoXX.cpp
	@$(MV)

obj/color.o: src/color.cpp src/color.h
	$(CC) -c src/color.cpp
	@$(MV)

obj/util.o: src/util.cpp src/util.h
	$(CC) -c src/util.cpp
	@$(MV)

obj/translator.o: src/translator.cpp src/translator.h src/gc.h
	$(CC) -c src/translator.cpp
	@$(MV)

obj/cleanliness.o: src/cleanliness.cpp src/cleanliness.h
	$(CC) -c src/cleanliness.cpp
	@$(MV)

obj/gzstream.o: src/gz/gzstream.C src/gz/gzstream.h	
	$(CFLAGS) -c src/gz/gzstream.C
	@$(MV)

# [Internal rules]

# Creating the object folder
obj/:
	@$(MK)

# Print info at compile start
start:
	@echo "";
	@echo "   ________________________________ \n";
	@echo "     <<< BUILDING SANS SERIF >>>  \n";
	@echo "   ________________________________";
	@echo "";

# Print info when done
done:
	@echo "";
	@echo "   _______________ \n";
	@echo "    <<< Done! >>> \n";
	@echo "   _______________";
	@echo "";

.PHONY: clean
clean:
	$(RM)




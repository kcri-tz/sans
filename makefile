# MAX. K-MER LENGTH, NUMBER OF FILES
CC = g++ -O3 -march=native -DmaxK=64 -DmaxN=64 -std=c++14

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
	MV = mv *.o obj/
	CP = cp makefile obj/makefile
	

endif

ifeq ("$(wildcard $(TD))", "")
    RM = @echo ""
endif

SANS:  start obj/ checkobs obj/main.o done
	$(CC) -o SANS obj/main.o obj/graph.o obj/kmer32.o obj/kmerXX.o obj/kmerAminoXX.o obj/kmerAmino12.o obj/color64.o obj/colorXX.o obj/util.o obj/translator.o obj/cleanliness.o obj/gzstream.o -lz $(BF)

obj/main.o: src/main.cpp src/main.h obj/translator.o obj/graph.o obj/util.o obj/cleanliness.o obj/gzstream.o
	$(CC) -c src/main.cpp
	$(MV)

obj/graph.o: src/graph.cpp src/graph.h obj/kmer32.o obj/kmerXX.o obj/kmerAmino12.o obj/kmerAminoXX.o obj/color64.o obj/colorXX.o
	$(CC) -c src/graph.cpp
	$(MV)

obj/kmer32.o: src/kmer32.cpp src/kmer32.h
	$(CC) -c src/kmer32.cpp
	$(MV)

obj/kmerXX.o: src/kmerXX.cpp src/kmerXX.h
	$(CC) -c src/kmerXX.cpp
	$(MV)

obj/kmerAmino12.o: src/kmerAmino12.cpp src/kmerAmino12.h obj/util.o
	$(CC) -c src/kmerAmino12.cpp
	$(MV)

obj/kmerAminoXX.o: src/kmerAminoXX.cpp src/kmerAminoXX.h obj/util.o
	$(CC) -c src/kmerAminoXX.cpp
	$(MV)

obj/color64.o: src/color64.cpp src/color64.h
	$(CC) -c src/color64.cpp
	$(MV)

obj/colorXX.o: src/colorXX.cpp src/colorXX.h
	$(CC) -c src/colorXX.cpp
	$(MV)

obj/util.o: src/util.cpp src/util.h
	$(CC) -c src/util.cpp
	$(MV)

obj/translator.o: src/translator.cpp src/translator.h src/gc.h
	$(CC) -c src/translator.cpp
	$(MV)

obj/cleanliness.o: src/cleanliness.cpp src/cleanliness.h
	$(CC) -c src/cleanliness.cpp
	$(MV)

obj/gzstream.o: src/gz/gzstream.C src/gz/gzstream.h	
	$(CFLAGS) -c src/gz/gzstream.C
	$(MV)

# [Internal rules]
obj/: makefile
	@$(MK)

# This rule checks if the makefile has changed
checkobs: makefile 
	@touch obj/makefile;
	@cmp -s "makefile" "obj/makefile"; \
	CONTINUE=$$?; \
	if [ $$CONTINUE -eq 0 ]; \
	then echo "\n   >>> NO CHANGES DETECTED: CONTINUING BUILD \n\n"; \
	else echo "\n   >>> MAKEFILE CHANGES DETECTED: REBUILDING \n\n" && $(RM) && $(MK) && $(CP); \
	fi


start: makefile
	@echo "";
	@echo "   ________________________________ \n";
	@echo "     <<< COMPILING SANS SERIF >>>  \n";
	@echo "   ________________________________";
	@echo "";

done: makefile
	@echo "";
	@echo "   _______________ \n";
	@echo "    <<< Done! >>> \n";
	@echo "   _______________";
	@echo "";

.PHONY: clean
clean:
	$(RM)




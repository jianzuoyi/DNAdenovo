all: kmerfreq
kmerfreq: main.cpp read.cpp kmer.cpp hashSet.cpp threads.cpp gzstream.cpp
	g++ -static -O3  -o $@ $^ -lz -lpthread 
clean:
	rm kmerfreq


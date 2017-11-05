all:  filter_data_parallel duplication
	
filter_data_parallel:  filter_low_quality.cpp filter_adapter.cpp filter_small_size.cpp filter_data_parallel.cpp gzstream.cpp
	g++  -O3 -static  -o $@  $^ -lz -lpthread

duplication: duplication.cpp gzstream.cpp
	g++ -O3 -static -o $@ $^ -lz 

clean:
	rm filter_data_parallel duplication

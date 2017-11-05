#ifndef __READ__H_
#define __READ__H_

#include <iostream>
#include <string>
#include <zlib.h>
#include <fstream>

#include "kmer.h"


extern uint64_t total_read;
extern uint64_t num_kmers;
extern uint64_t used_base;
extern uint64_t num_reads;

extern int Read_len;
extern int Ahead;
extern int Delete;
extern uint64_t genome;

void LoadList(string &files, bit8_t* &freq, HashSet * &kmerHash);
//void LoadList(string &files, bit8_t* &freq);
void LoadList(string &files);
#endif


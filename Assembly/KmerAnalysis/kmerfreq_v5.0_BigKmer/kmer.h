#ifndef __KMER__H_
#define __KMER__H_

#include <fstream>
#include <zlib.h>
#include <iostream>

#include "hashSet.h"
#include "gzstream.h"
//#include "global.h"

#define UPLIMIT 255
typedef unsigned char bit8_t;

extern int Kmer;
extern int Min_number;
extern int times_analysis;
extern int KPS;
extern int prefix_bit_num;
extern uint64_t num_nodes;
extern uint64_t num_kmers;


//�ɣ��ãǣԵ�ASCII�뵽�������������Զ������Сд
//256����ĸ��alphabet����,��8bit��char�ʹ洢,A=a=0,C=c=1,G=g=2,T=t=3,��������ĸ��Ϊ4
const char alphabet[256] =
{
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 
 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};

const char sbin[4] ={
	'A', 'C', 'G', 'T'
};

typedef struct KmerElement{
        char kmerSeq[35];
        uint32_t kmerFreq;
};

void Getfreqy(bit8_t *&freq, uint64_t total, int* &freqy);
int GetDepth(uint64_t* &freqy, int &min_depth, int &min_freq);

void Read2Kmer(int buffId, uint64_t readId, int threadId);
string Kmer2Seq(uint64_t mer, int kmer);
void CombineKmer(bit8_t* &freq, uint64_t total, HashSet * &kmerHash);
void CombineKmer(bit8_t* &freq, uint64_t total, HashSet** &kmerHashSet);
void CombineKmer(bit8_t* &freq, uint64_t total);

void HandleKmer0(int buffId, uint64_t index, int threadId);
void HandleKmer1(int buffId, uint64_t index, int threadId);

void WriteKmer(bit8_t* &freq, uint64_t total);
void WriteKmer(HashSet* &kmerHash); //out in *.depth.gz
void WriteKmer1(HashSet* &kmerHash); //out in *.depth
void WriteKmer1(HashSet** &kmerHashSet);
void WriteKmer(bit8_t* &freq, uint64_t total, HashSet * &kmerHash);
void WriteKmer(bit8_t* &freq, uint64_t total, HashSet** &kmerHashSet);

void ReadKmer();
void ReadKmer(string &file, uint8_t *newfreq);

void outputKmer(HashSet**& kmerHashSet, int threadNum, uint64_t*& freqy);

#endif

#include <iostream>
#include <string>
#include <sys/time.h>
#include <fstream>
#include <vector>
#include <algorithm>


#include "read.h"
#include "kmer.h"
#include "threads.h"

using namespace std;

int Kmer = 17; //14mer
int Min_number = 1; //设置输出的最小出现个数
int Delete=0;  //设置末尾删除的长度
int Read_len = 1; //设置reads的长度，小于该长度的忽略，大于的截短
int Ahead=0;
uint64_t genome=0;
int mode=0; //simple mode 0, full mode 1.
int output_type=0; //1 is output the frequence high depth(>255)a file and low depth a file. 2 is output them together. 0 is not output.

bit8_t* freq = NULL;

int threadNum = 1;		//thread number
int maxReadLen = 150;		//maximum read length
int bufferNum = 2;		//number of buffer
int bufferSize = 5000000;	//maximum kmer number in each buffer
uint64_t** kmerBuffer = NULL;	//store kmers
uint64_t* kmerNum = NULL;	//current kmer number in each buffer
uint64_t** hcBuffer = NULL;	//store hashcode of kmers
uint64_t** indexArr = NULL;	//store substript of first kmer in each read
uint64_t maxReadNum = 0;	//maximu read number in each buffer
uint64_t* readNum = NULL;	//current read number in each buffer
string** readBuffer = NULL;	//store reads
HashSet** kmerHashSet = NULL;	//kmer hashs for all threads
int** isBuffFull = NULL;	//full state of buffers
int** isBuffEmpty = NULL;	//empty of buffers
int* threadSIG = NULL;		//signals of threads
int* pthreadId = NULL;		//id of thread (eg. 0~7 when threadNum=8)
pthread_t* threads = NULL;	

double load_factor=0.75;	//load factor of hash table
int initSize=1024*1024;		//initialized hash table size

string file_list;

uint64_t num_kmers=0; //kmer number
uint64_t num_reads=0; //real used read number
uint64_t total_read=0; //total input read number
uint64_t num_nodes=0; //unique kmer number
uint64_t used_base=0; //used base number

void usage(void)
{
	cout<<"\nUsage:  kmerfreq [options]\n"
	    <<"version: 1.0 2010-5-10 \n"
	    <<"Contact: liubinghang@genomics.org.cn\n"
		<<"	-k <int>	kmer size(13~17, need 4^k bytes memory), default="<<Kmer<<"\n"
		<<"	-r <int>	read length used to get kmers, default=read's real length\n"
		<<"	-a <int>	ignored length of the beginning of a read, default="<<Ahead<<"\n"
		<<"	-n <int>	minimal appeared times, default="<<Min_number<<"\n"
		<<"	-d <int>	ignored length of the end of a read, default="<<Delete<<"\n"
		<<"	-g <int>	total bases used to get kmers, default=all input bases\n"
		<<"	-l <string>	read file list\n"
		<<"	-i <int>	initial size of hash table, default="<<initSize<<"\n"
		<<"	-m <int>	simple mode(0, depth larger than 255 will be treated as 255), full mode(1, record real depth) "<<mode<<"\n"
		<<"	-o <int>	output type(1: output low and high depth(>=255) in two file, 2: output them together, 0: do not output depth file), default= "<<output_type<<"\n"
		<<"	-t <int>	thread number, default=" << threadNum << "\n"
		<<"	-L <int>	maximum read length, default=" << maxReadLen << "\n"
		<<"	-h		output help information to screen\n\n"
		<<"Example: \nkmerfreq -k 17 -l fq.list >17mer.freq 2>17mer.log;\n"
		<<"kmerfreq -k 17 -m 1 -o 1 -l fq.list >17mer.freq 2>17mer.log;\n\n"
		<<"Attension: Please don't set -d and -r at the same time.\n\n";
	exit(0);
}

int main(int argc, char *argv[])
{

	if (argc < 2) usage();

	int c;
	while((c=getopt(argc, argv, "k:r:n:a:d:g:l:i:m:o:t:L:h")) !=-1) {
		switch(c) {
			case 'k': Kmer=atoi(optarg); break;
			case 'r': Read_len=atoi(optarg); break;
			case 'n': Min_number=atoi(optarg); break;
			case 'a': Ahead=atoi(optarg);break;
			case 'd': Delete=atoi(optarg);break;
			case 'g': genome=strtol(optarg, NULL, 10);break;
			case 'l': file_list = optarg; break;
			case 'i': initSize = atoi(optarg);break;
			case 'm': mode = atoi(optarg);break;
			case 'o': output_type = atoi(optarg); break;
			case 't': threadNum = atoi(optarg); break;
			case 'L': maxReadLen = atoi(optarg); break;
			case 'h': usage(); break;
			default: usage();
		}
	}

	if (Kmer < 13 || Kmer > 17)
	{
		cerr << "Kmer should be between 13 and 17!\n";
		exit(3);
	}

	if ((Read_len > 1) && (Read_len < Kmer))
	{
		cerr << "r(" << Read_len << ") is smaller than Kmer(" << Kmer << "). Please check the parameter!\n";
		exit(3);
	}

	if (threadNum < 1)
	{
		threadNum = 1;
	}

	if (maxReadLen <= 30)
	{
		maxReadLen = 100;
	}

	if (Read_len > 1)
	{
		maxReadLen = Read_len;
	}

	kmerBuffer = new uint64_t*[bufferNum];
	kmerNum = new uint64_t[bufferNum];
	hcBuffer = new uint64_t*[bufferNum];
	indexArr = new uint64_t*[bufferNum];
	readBuffer = new string*[bufferNum];
	readNum = new uint64_t[bufferNum];
	isBuffFull = new int*[bufferNum];
	isBuffEmpty = new int*[bufferNum];

	maxReadNum = uint64_t(bufferSize/(maxReadLen-Kmer+1));

	for (int i=0; i<bufferNum; i++)
	{
		kmerNum[i] = 0;
		readNum[i] = 0;

		kmerBuffer[i] = new uint64_t[bufferSize];
		hcBuffer[i] = new uint64_t[bufferSize];

		indexArr[i] = new uint64_t[maxReadNum];
		readBuffer[i] = new string[maxReadNum];

		isBuffFull[i] = new int[threadNum];
		isBuffEmpty[i] = new int[threadNum];
		for (int j=0; j<threadNum; j++)
		{
			isBuffFull[i][j] = 0;
			isBuffEmpty[i][j] = 1;
		}
	}
	
	threadSIG = new int[threadNum];
	pthreadId = new int[threadNum];
	threads = new pthread_t[threadNum];
	kmerHashSet = new HashSet*[threadNum];
	for (int i=0; i<threadNum; i++)
	{
		threadSIG[i] = 0;
		pthreadId[i] = i;
		int it;
		if ((it=pthread_create(threads+i, NULL, threadRoutine, (void*)(pthreadId+i))) != 0)
		{
			cerr << "Failed to create thread!\n";
			exit(4);
		}

		if (mode == 1)
		{
			kmerHashSet[i] = init_hashset(initSize, load_factor);
		}
	}

	cerr << threadNum << " threads created!\n";

	time_t start_time = time(NULL);
	uint64_t total=0;  //理论上最大的kmer个数 

	int *freqy = NULL; //frequency 

	for(int i=0; i<Kmer; i++) {
		total=(total<<2)|0x3;
	}
	freq=new bit8_t[total]; //动态内存申请

	for(uint64_t i=0; i<total; i++) {
		freq[i]=0; 
	}
	freqy=new int[256];
	for(int i=1;i<256;i++){
		freqy[i]=0;
	}

	LoadList(file_list, freq);

	for (int i=0; i<threadNum; i++)
	{
		threadSIG[i] = 1;
	}

	cerr << "Waitint for threads to exit...\n";
	for (int i=0; i<threadNum; i++)
	{
		pthread_join(threads[i], NULL);
	}

	cerr << "All threads finished!\nUsed time: " << double(time(NULL) - start_time)/60.0 << " minutes.\n";


	if (mode == 0)
	{
		CombineKmer(freq, total);
		if(output_type != 0)
		{
			WriteKmer(freq, total);
		}
	}
	else
	{
		CombineKmer(freq, total, kmerHashSet);
		if(output_type ==1)
		{
			WriteKmer(freq, total);
			WriteKmer1(kmerHashSet);
		}
		else if(output_type ==2)
		{
			WriteKmer(freq, total, kmerHashSet);
		}
	}

	Getfreqy(freq, total, freqy);

	int Min_depth, Min_freq;
	int Expect_depth=GetDepth(freqy, Min_depth, Min_freq);

	//output the result table.
	cerr<<"kmer\tkmer_num\tpkdepth\tgenome_size\tused_base\tused_read\tX\taverage_read_len\tnode_num\n"
	    <<Kmer<<"\t"<<num_kmers<<"\t"<<Expect_depth<<"\t"<<num_kmers/Expect_depth<<"\t"
	    <<used_base<<"\t"<<num_reads<<"\t"<<1.00*used_base/num_kmers*Expect_depth<<"\t"
	    <<used_base/total_read<<"\t"<<num_nodes<<endl;

	//remove the low depth part
	for(int i=1;i< Min_depth;i++){
		num_kmers=num_kmers-freqy[i]+Min_freq;
	}
	cerr<<"modified kmer_num is: "<<num_kmers<<"\t"<<num_kmers/Expect_depth<<endl;

	for (int i=0; i<bufferNum; i++)
	{
		delete[] kmerBuffer[i];
		delete[] hcBuffer[i];
		delete[] indexArr[i];
		delete[] readBuffer[i];
		delete[] isBuffFull[i];
		delete[] isBuffEmpty[i];
	}

	if (mode == 1)
	{
		for (int i=0; i<threadNum; i++)
		{
			free_hash(kmerHashSet[i]);
		}
	}

	delete[] kmerHashSet;
	delete[] kmerBuffer;
	delete[] hcBuffer;
	delete[] indexArr;
	delete[] readBuffer;
	delete[] threadSIG;
	delete[] pthreadId;
	delete[] threads;
	delete[] isBuffFull;
	delete[] isBuffEmpty;
	delete[] readNum;
	delete[] kmerNum;

	delete []freqy;
	delete []freq;
	time_t end_time = time(NULL);
	cerr<<"\nAll done! Using time "<<double(end_time-start_time)/60.0<<"mins.\nThank you!\n"<<endl;

	return 0;
}

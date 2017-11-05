#include "kmer.h"

extern int threadNum;
extern HashSet** kmerHashSet;
extern uint64_t** indexArr;
extern string** readBuffer;
extern uint64_t** hcBuffer;
extern uint64_t** kmerBuffer;
extern bit8_t* freq;
extern int maxReadLen;


//convert sequence to kmer-bit，64bit最多能装32mer
//需要事先确定序列中只含有ACGT(或acgt)4个碱基
uint64_t Seq2Kmer(string &seq)
{
	uint64_t mer=0;
	for(int i=0; i<seq.size(); i++) {
		mer=(mer<<2)|alphabet[seq[i]];
	}
	return mer;
}

//convert kmer-bit to kmer-sequence
//64bit最多能装32mer,此处必须给定kmer
string Kmer2Seq(uint64_t mer, int kmer)
{
	string seq;
	for(int i=0; i<kmer; i++) {
		seq.push_back(sbin[(mer>>(kmer-1-i)*2)&0x3]);
	}
	return seq;
}

uint64_t get_rev_seq(uint64_t seq, int seq_size){
        seq = ~seq;
        //cout << seq << "\n";
        seq = ((seq & 0x3333333333333333LLU)<<  2) | ((seq & 0xCCCCCCCCCCCCCCCCLLU)>>  2);
        seq = ((seq & 0x0F0F0F0F0F0F0F0FLLU)<<  4) | ((seq & 0xF0F0F0F0F0F0F0F0LLU)>>  4);
        seq = ((seq & 0x00FF00FF00FF00FFLLU)<<  8) | ((seq & 0xFF00FF00FF00FF00LLU)>>  8);
        seq = ((seq & 0x0000FFFF0000FFFFLLU)<< 16) | ((seq & 0xFFFF0000FFFF0000LLU)>> 16);
        seq = ((seq & 0x00000000FFFFFFFFLLU)<< 32) | ((seq & 0xFFFFFFFF00000000LLU)>> 32);
        return seq >> (64 - (seq_size<<1));
}


string reverse_complement (string &seq)
{	
	string nseq;
	for (int i=seq.size()-1; i>=0; i--)
	{	char base;
		if (seq[i] == 'A')
			base = 'T';
		if (seq[i] == 'C')
			base = 'G';
		if (seq[i] == 'G')
			base = 'C';
		if (seq[i] == 'T')
			base = 'A';
		nseq.push_back(base);
	}
	return nseq;
}

void Read2Kmer(int buffId, uint64_t readId, int threadId)
{
	uint64_t index = indexArr[buffId][readId];

	for (string::size_type pos=0; pos<readBuffer[buffId][readId].size()-Kmer+1; pos++)
	{
		string kmerSeq = readBuffer[buffId][readId].substr(pos, Kmer);
		uint64_t kmerBit = Seq2Kmer(kmerSeq);

		uint64_t hashcode = HashCode(kmerBit);
		hcBuffer[buffId][index] = hashcode;
		kmerBuffer[buffId][index++] = kmerBit;
	}
}

uint64_t Minbit(uint64_t &mer, uint64_t &rcmer)
{
	return mer > rcmer ? rcmer : mer;
}


void CombineKmer(bit8_t* &freq, uint64_t total)
{
	for(uint64_t i=0; i<total; i++)
	{
		if(freq[i] > 0)
		{
			num_nodes++;
			uint64_t rcmer=get_rev_seq( i, Kmer);
			if (i == rcmer)
			{
				continue;
			}

			if(freq[i] + freq[rcmer] < UPLIMIT)
			{
				if(i > rcmer)
				{
					freq[rcmer]+=freq[i];
					freq[i]=0;
				}else{
					freq[i]+=freq[rcmer];
					freq[rcmer]=0;
				}
			}else{
				if(i > rcmer)
				{
					freq[rcmer]=UPLIMIT;
					freq[i]=0;
				}else{
					freq[i]=UPLIMIT;
					freq[rcmer]=0;
				}
			}
		}
	}
}

void CombineKmer(bit8_t* &freq, uint64_t total, HashSet** &kmerHashSet)
{
	for(uint64_t i=0; i<total; i++) 
	{
		if((int)freq[i] > 0) 
		{
			num_nodes++;

			uint64_t hashcode = HashCode(i);
			int threadId = hashcode%threadNum;
			uint64_t rcmer = get_rev_seq( i, Kmer);
			if (i == rcmer)
			{
				continue;
			}
			uint64_t rcHashcode = HashCode(rcmer);
			int rcThreadId = rcHashcode%threadNum;

			if ((int)freq[rcmer] >  0)
			{	
				if ((int)freq[i] + (int)freq[rcmer] < UPLIMIT)
				{	
					if(i > rcmer)
					{
						freq[rcmer] += freq[i];
						freq[i]=0;
					}else{
						freq[i] += freq[rcmer];
						freq[rcmer]=0;
					}

				}else{
					Entity* entity = new Entity;
					entity->kmer = i;
					uint64_t idx = get_hashset(kmerHashSet[threadId], entity);
					if(idx == kmerHashSet[threadId]->size) //kmer not exits
					{
						entity->kmer = rcmer;
						idx = get_hashset(kmerHashSet[rcThreadId], entity);

						if(idx == kmerHashSet[rcThreadId]->size) //both kmer and rckmer do not exist in hashset. add one.
						{
							entity->freq = (uint64_t)freq[i] + (uint64_t)freq[rcmer];
							if(i > rcmer)
							{
								entity->kmer = rcmer;
								freq[rcmer] = UPLIMIT;
								freq[i]=0;
								add_hashset(kmerHashSet[rcThreadId], entity);
							}else{
								entity->kmer = i;
								freq[i]=UPLIMIT;
								freq[rcmer]=0;
								add_hashset(kmerHashSet[threadId], entity);
							}
							delete entity;

						}else{
							kmerHashSet[rcThreadId]->array[ idx ].freq += (uint64_t)freq[i];
							freq[i]=0;

							if(i < rcmer)
							{
								entity->kmer = i;
								entity->freq = kmerHashSet[rcThreadId]->array[ idx ].freq;
								add_hashset(kmerHashSet[threadId], entity);
								entity->kmer = rcmer;
								int index=delete_hashset(kmerHashSet[rcThreadId], entity);
								freq[i] =UPLIMIT;
								freq[rcmer] =0;
							}

							delete entity;
						}

					}else{ //kmer exits
						entity->kmer = rcmer;
						uint64_t idx2 = get_hashset(kmerHashSet[rcThreadId], entity);

						if (idx2 == kmerHashSet[rcThreadId]->size)
						{

							if (i > rcmer)
							{
								entity->freq = kmerHashSet[threadId]->array [idx].freq + (uint64_t)freq[rcmer];
								add_hashset(kmerHashSet[rcThreadId], entity);
								freq[rcmer] = UPLIMIT;

								entity->kmer = i;
								delete_hashset(kmerHashSet[threadId], entity);
								freq[i] = 0;
							}
							else
							{
								kmerHashSet[threadId]->array [idx].freq += (uint64_t)freq[rcmer];
								freq[i] = UPLIMIT;
								freq[rcmer] = 0;
							}
						}
						else
						{
							if (i > rcmer)
							{
								kmerHashSet[rcThreadId]->array[idx2].freq += kmerHashSet[threadId]->array[idx].freq;
								freq[rcmer] = UPLIMIT;

								entity->kmer = i;
								delete_hashset(kmerHashSet[threadId], entity);
								freq[i] = 0;
							}
							else
							{
								kmerHashSet[threadId]->array [idx].freq += kmerHashSet[rcThreadId]->array[idx2].freq;
								freq[i] = UPLIMIT;

								entity->kmer = rcmer;
								delete_hashset(kmerHashSet[rcThreadId], entity);
								freq[rcmer] = 0;
							}
						}
						delete entity;
					} 
				}//end of more than 255. 

			}//reverse less than 1

		} //i freq less than 1
	}//for end
}

void WriteKmer(bit8_t* &freq, uint64_t total)
{
	uint64_t bufUnit = 2000000;
	gzFile wzip;
	wzip=gzopen("kmer.low.depth.gz", "wb");
		
	for (uint64_t i=0; i<total; i+=bufUnit)
	{	uint64_t bufSize = (total - i > bufUnit) ? bufUnit : (total - i);
		gzwrite(wzip, (voidpc)(freq+i), bufSize);
	}
	gzclose(wzip);
	cerr<<"Low Depth output finished!\n";
}


void WriteKmer1(HashSet** &kmerHashSet)
{
	ofstream outFile("kmer.high.depth");
	
	cerr<<"Output the kmer frequence without compression\n";

	for (int t=0; t<threadNum; t++)
	{
		for (uint64_t i=0; i<kmerHashSet[t]->size; i++)
		{
			if ((is_entity_null(kmerHashSet[t]->nul_flag, i) == 0) && (is_entity_delete(kmerHashSet[t]->del_flag, i) == 0))
			{
				outFile<<Kmer2Seq( kmerHashSet[t]->array[i].kmer, Kmer)<<"\t"<<int (kmerHashSet[t]->array[i].freq)<<endl;
			}
		}
	}
	outFile.close();

}

void WriteKmer(bit8_t* &freq, uint64_t total, HashSet** &kmerHashSet)
{
	gzFile wzip;
	wzip=gzopen("kmer.total.depth.gz", "wb");
	int buffSize= 5; //five number once
	int *buff;
	buff=new int[buffSize];
	for(int j=0; j< buffSize; j++)
		buff[j]=0;

	int now_num=0;
	uint64_t i=0;
        for(; i< total; i++)
        {
		if(i - now_num == buffSize)
		{
			gzwrite(wzip, buff, buffSize);
			now_num= i;
		}
                if(freq[i] < UPLIMIT)
                {
                        cerr<<Kmer2Seq(i, Kmer)<<"\t"<<int (freq[i])<<endl;
			buff[i - now_num]=int (freq[i]);
                }else{
                        Entity* entity= new Entity;
                        entity->kmer=i;
			uint64_t hashcode = HashCode(i);
			int threadId = hashcode%threadNum;
                        uint64_t idx = get_hashset(kmerHashSet[threadId], entity);
                        if(idx != kmerHashSet[threadId]->size)
                        {
                                cerr<<Kmer2Seq(i, Kmer)<<"\t"<<int (freq[i])<<"\t"<<int(kmerHashSet[threadId]->array[idx].freq)<<endl;
				buff[i - now_num] = int(kmerHashSet[threadId]->array[idx].freq);
                        }else{
                                cerr<<"ERROR for frequence\n";
				buff[i - now_num] = int (freq[i]);
                        }
                        delete entity;
                }
        }
	int dis=i - now_num;
	if(dis >= 0)
	{
		int *subbuff=new int [dis +1];
		for(int j=0; j< dis+1; j++)
			subbuff[j]=buff[j];
		gzwrite(wzip, subbuff, dis);
		delete []subbuff;
	}
	delete []buff;
        cerr<<"All Kmer write done!\n";
        gzclose(wzip);
}

void ReadKmer()
{
	gzFile zip;
	zip=gzopen("kmer.freq.gz", "rb");
	int c;
	cout<<"Start read file "<<endl;
	while((c=gzgetc(zip)) != EOF)
	{
		cout<<int (c)<<" ";
	}
	cout<<endl;
	cerr<<"All kmer read done!\n";
	gzclose(zip);
}

void ReadKmer(string &file, uint8_t *newfreq) //make sure to initial newfreq before useing this function
{
	uint64_t bufUnit = 2000000;
	gzFile zip;
	zip=gzopen(file.c_str(), "rb");
	uint64_t i=0;
	uint64_t count_num = 0;
	do
	{	
		count_num = gzread(zip, (voidp)(newfreq+i), bufUnit);
		i += bufUnit;
		
	}
	while (count_num == bufUnit);
	cerr << "gz decompress generated" << endl;
	gzclose(zip);
}


void Getfreqy(bit8_t *&freq, uint64_t total, int* &freqy)
{
	for(uint64_t i=0; i<=total; i++) {
		if(freq[i] >= Min_number) {
			freqy[int(freq[i])]++;
		}
	}
}

int GetDepth(int* &freqy, int &min_depth, int &min_freq)
{
	//find the peek position and output the frequence file.
	int max_freq=0;
	int max_depth=1; 
	int find_max=0; int find_min=0;

	for(int i=1;i<256;i++){
		cout<<i<<"\t"<<freqy[i]<<"\t"<<freqy[i]*100.00/num_nodes<<"\t"<<i*freqy[i]<<"\n"; 
		if(freqy[i]<freqy[i+1] && find_min==1){max_freq=freqy[i+1];max_depth=i+1;find_max=1;}
		if(freqy[i]>freqy[i+1] && find_max==0){min_freq=freqy[i+1];min_depth=i+1;find_min=1;}
		if(freqy[i]>freqy[i+1] && find_max==1){find_min=0;}	
	}
	
	//check the peek
	if(freqy[max_depth+1]*(max_depth+1) >= freqy[max_depth+2]*(max_depth+2) && freqy[max_depth+1]*(max_depth+1) >= freqy[max_depth]*max_depth ){
		cerr<<"The peek position is true"<<endl;
	}else{
		cerr<<"Please check the peek position carefully"<<endl;
	}
	return max_depth;
}

void HandleKmer0(int buffId, uint64_t index, int threadId)
{
	uint64_t kmerBit = kmerBuffer[buffId][index];

	if (freq[kmerBit] < UPLIMIT)
	{
		freq[kmerBit]++;
	}
}

void HandleKmer1(int buffId, uint64_t index, int threadId)
{
	uint64_t kmerBit = kmerBuffer[buffId][index];
	if (freq[kmerBit] < UPLIMIT)
	{
		freq[kmerBit]++;
	}
	else
	{
		Entity* entity = new Entity;
		entity->kmer = kmerBit;
		entity->freq = 256;
		uint64_t idx = get_hashset(kmerHashSet[threadId], entity);
		if(idx == kmerHashSet[threadId]->size)
			add_hashset(kmerHashSet[threadId], entity);
		else
			kmerHashSet[threadId]->array[idx].freq++;
		delete entity;
		
	}
}

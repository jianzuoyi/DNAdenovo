/****************************************************
*	author: Li Zhenyu (lizhenyu@genomics.cn)
*		Fan Wei (fanw@genomics.cn)
*	version: 3.0
*
****************************************************/

#include"hashSet.h"

int Kmer = 50;
uint64_t init_size = 1024*1024;
float load_factor = 0.75;

time_t start_time;
ofstream log_outfile, kmer_outfile, kmerFreq_outfile, unique_outfile, repeat_outfile;
uint32_t useq_count=0, rseq_count=0; //count uinque sequence number and repeat sequence number
uint64_t max_unique_len=0, max_repeat_len=0; // maximum uinque sequence length, maximum repeat sequence length
string::size_type useq_len=0, rseq_len=0; //total unique sequence length, total repeat sequence length

char alphabet[256] =
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

//convert sequence to kmer-number
void Seq2Kmer(uint64_t* kmer, string& seq)
{
	if (seq.size() > 32)
	{
		int i=0;
		for ( ; i<seq.size()-32; i++)
		{
			kmer[1] = (kmer[1]<<2)|alphabet[seq[i]];
		}

		for ( ; i<seq.size(); i++)
		{
			kmer[0] = (kmer[0]<<2)|alphabet[seq[i]];
		}
	}
	else
	{
		for (int i=0; i<seq.size(); i++)
		{
			kmer[0] = (kmer[0]<<2)|alphabet[seq[i]];
		}
	}
}

uint64_t get_rev_seq(uint64_t seq, int seq_size)
{
	seq = ~seq;
	seq = ((seq & 0x3333333333333333LLU)<<  2) | ((seq & 0xCCCCCCCCCCCCCCCCLLU)>>  2);
	seq = ((seq & 0x0F0F0F0F0F0F0F0FLLU)<<  4) | ((seq & 0xF0F0F0F0F0F0F0F0LLU)>>  4);
	seq = ((seq & 0x00FF00FF00FF00FFLLU)<<  8) | ((seq & 0xFF00FF00FF00FF00LLU)>>  8);
	seq = ((seq & 0x0000FFFF0000FFFFLLU)<< 16) | ((seq & 0xFFFF0000FFFF0000LLU)>> 16);
	seq = ((seq & 0x00000000FFFFFFFFLLU)<< 32) | ((seq & 0xFFFFFFFF00000000LLU)>> 32);
	return seq >> (64 - (seq_size<<1));
}

void getRCKmer(uint64_t *RCKmer, uint64_t *kmer, int Kmer)
{
	if (Kmer <= 32)
	{
		RCKmer[1] = 0;
		RCKmer[0] = get_rev_seq(kmer[0], Kmer);
	}
	else
	{
		RCKmer[1] = get_rev_seq(kmer[0], 32);
		RCKmer[0] = get_rev_seq(kmer[1], Kmer-32);

		RCKmer[0] |= RCKmer[1] << ((Kmer-32)<<1);
		RCKmer[1] = RCKmer[1] >> ((64-Kmer)<<1);
	}
}

bool isLargerKmer(uint64_t *kmer, uint64_t *RCKmer)
{
	if (kmer[1] > RCKmer[1])
	{
		return true;
	}
	else if (kmer[1] < RCKmer[1])
	{
		return false;
	}
	else
	{
		return (kmer[0] > RCKmer[0]);
	}
}

inline void Usage()
{
	cerr << "Usage: /path/kmer_scan [option] *.fa\n"
		<< "\t-k <int>		Kmer length(<=60), default=" << Kmer << "\n"
		<< "\t-m <int>		initial size of hash table, default=" << init_size << "\n"
		<< "\t-l <float>	load factor of hash table, default=" << load_factor << "\n"
		<< "\t-h		output this information\n"
		<< "example:\n\t/path/kmer_scan ./test.fa\n\t/path/kmer_scan -k 55 -m 16700000 ./test.fa\n";
}

//When reaching the end of unique/repeat region, update related variants and output them
inline void reach_end(string& seq, string::size_type& start_pos, string::size_type& end_pos, string& id,
		uint32_t& count, uint64_t& max_len, string::size_type& length, ofstream& outfile)
{
	string subseq = seq.substr(start_pos, end_pos-start_pos);
	max_len = (max_len<(end_pos-start_pos)) ? (end_pos-start_pos) : max_len;
	length += end_pos - start_pos;
	outfile << ">" << count << "\t" << start_pos << "\t" << end_pos-1
			<< "\t" << subseq.size() << "\t" << id << "\n" << subseq << "\n";
}


//fill hash table
void fill_hash_table(HashSet *Hash, string& seq)
{
	//check if there is 'N' in sequence
	if (seq.find('N') != string::npos) { log_outfile << "'N' is not allowed in sequence\n"; exit(1); }

	string kmer_seq;
	uint64_t leftBase = 4;
	uint64_t rightBase = 4;
	uint64_t temp;

	if (seq.size() > Kmer)
	{
		rightBase = alphabet[seq[Kmer]];
	}

	uint64_t kmer[2] = {0, 0};
	uint64_t RCKmer[2] = {0, 0};

	kmer_seq = seq.substr(0, Kmer);

	Seq2Kmer(kmer, kmer_seq);
	getRCKmer(RCKmer, kmer, Kmer);

	if (isLargerKmer(kmer, RCKmer))
	{
		kmer[0] = RCKmer[0];
		kmer[1] = RCKmer[1];

		temp = leftBase;
		leftBase = rightBase;
		rightBase = 3 - temp;
	}

	Entity entity;
	entity.high = kmer[1];
	entity.low  = kmer[0];
	entity.freq = 1;
	entity.leftBranch = 0;
	entity.rightBranch = 0;
	entity.leftBase = leftBase;
	entity.rightBase = rightBase;

	add_hashset(Hash, &entity);

	
	string::size_type pos;
	for (pos=1; pos<(seq.size()-Kmer); pos++)
	{
		leftBase = alphabet[seq[pos-1]];
		rightBase = alphabet[seq[pos+Kmer]];

		kmer_seq = seq.substr(pos, Kmer);

		kmer[0] = 0;
		kmer[1] = 0;
		RCKmer[0] = 0;
		RCKmer[1] = 0;

		Seq2Kmer(kmer, kmer_seq);
		getRCKmer(RCKmer, kmer, Kmer);

		if (isLargerKmer(kmer, RCKmer))
		{
			kmer[0] = RCKmer[0];
			kmer[1] = RCKmer[1];

			temp = leftBase;
			leftBase = 3 - rightBase;
			rightBase = 3 - temp;
		}

		entity.high = kmer[1];
		entity.low  = kmer[0];
		entity.freq = 1;
		entity.leftBranch = 0;
		entity.rightBranch = 0;
		entity.leftBase = leftBase;
		entity.rightBase = rightBase;

		add_hashset(Hash, &entity);
	}

	if (seq.size() >= pos+Kmer)
	{
		leftBase = alphabet[seq[pos-1]];
		rightBase = 4;

		kmer_seq = seq.substr(pos, Kmer);

		kmer[0] = 0;
		kmer[1] = 0;
		RCKmer[0] = 0;
		RCKmer[1] = 0;

		Seq2Kmer(kmer, kmer_seq);
		getRCKmer(RCKmer, kmer, Kmer);

		if (isLargerKmer(kmer, RCKmer))
		{
			kmer[0] = RCKmer[0];
			kmer[1] = RCKmer[1];

			temp = leftBase;
			leftBase = rightBase;
			rightBase = 3 - temp;
		}

		entity.high = kmer[1];
		entity.low  = kmer[0];
		entity.freq = 1;
		entity.leftBranch = 0;
		entity.rightBranch = 0;
		entity.leftBase = leftBase;
		entity.rightBase = rightBase;

		add_hashset(Hash, &entity);
	}
}
	
	
//get unique and repeat sequences and output them
void fish_sequence(HashSet *Hash, string &seq, string &id)
{
	string::size_type pos=0, start_pos, end_pos;
	uint64_t freq=0;
	int repeat2unique_flag = 0; //record the switch from repeat sequence to unique sequence

	uint64_t leftBranch = 0;
	uint64_t rightBranch = 0;

	string kmer_seq = seq.substr(pos, Kmer);
	uint64_t kmer[2] = {0, 0};
	uint64_t RCKmer[2] = {0, 0};
	Seq2Kmer(kmer, kmer_seq);
	getRCKmer(RCKmer, kmer, Kmer);

	if (isLargerKmer(kmer, RCKmer))
	{
		kmer[0] = RCKmer[0];
		kmer[1] = RCKmer[1];
	}

	freq = get_freq(Hash, kmer, &leftBranch, &rightBranch);

	while (pos < (seq.size()-Kmer+1))
	{
		if (freq == 1)
		{
			repeat2unique_flag = 0;

			start_pos = pos;
			end_pos = start_pos + Kmer;
			useq_count++;

			while (1)
			{
				pos++;

				if (pos == seq.size()-Kmer+1) 
				{/* reached the end of genome sequence */

					if (start_pos > 0)
					{
						start_pos--;
					}

					reach_end(seq, start_pos, end_pos, id, useq_count, 
							max_unique_len, useq_len, unique_outfile);
					break; 
				}

				kmer_seq = seq.substr(pos, Kmer);
				kmer[0] = 0; 
				kmer[1] = 0;
				Seq2Kmer(kmer, kmer_seq);
				RCKmer[0] = 0;
				RCKmer[1] = 0;
				getRCKmer(RCKmer, kmer, Kmer);

				if (isLargerKmer(kmer, RCKmer))
				{
					kmer[0] = RCKmer[0];
					kmer[1] = RCKmer[1];
				}

				freq = get_freq(Hash, kmer, &leftBranch, &rightBranch);

				if (freq == 1)
				{/* still unique kmer */
					end_pos++;
				}
				else if (freq > 1)
				{/* met repeat kmer */
					if (start_pos > 0)
					{
						start_pos--;
					}

					end_pos++;
					reach_end(seq, start_pos, end_pos, id, useq_count, 
							max_unique_len, useq_len, unique_outfile);
					break;
				}
			}
		}

		if (freq > 1)
		{
			repeat2unique_flag = 0;

			start_pos = pos;
			end_pos = start_pos + Kmer;

			while (1)
			{
				pos++;

				if (pos == seq.size()-Kmer+1) 
				{/* reached the end of genome sequence */

					rseq_count++;

					reach_end(seq, start_pos, end_pos, id, rseq_count, 
							max_repeat_len, rseq_len, repeat_outfile);
					break;
				}

				if ((rightBranch == 1 || leftBranch == 1) && (pos > start_pos+1))
				{
					rseq_count++;

					reach_end(seq, start_pos, end_pos, id, rseq_count, max_repeat_len, rseq_len, repeat_outfile);
					start_pos = pos - 1;
					end_pos = start_pos + Kmer;
				}

				kmer_seq = seq.substr(pos, Kmer);
				kmer[0] = 0;
				kmer[1] = 0;
				Seq2Kmer(kmer, kmer_seq);

				RCKmer[0] = 0;
				RCKmer[1] = 0;
				getRCKmer(RCKmer, kmer, Kmer);

				if (isLargerKmer(kmer, RCKmer))
				{
					kmer[0] = RCKmer[0];
					kmer[1] = RCKmer[1];
				}

				freq = get_freq(Hash, kmer, &leftBranch, &rightBranch);


				if (freq > 1)
				{/* still repeat kmer */
					end_pos++;
				}
				else if (freq == 1)
				{/* met unique kmer */

					if (pos > start_pos + 1)
					{
						reach_end(seq, start_pos, end_pos, id, rseq_count, 
								max_repeat_len, rseq_len,repeat_outfile);
					}
					repeat2unique_flag = 1;

					break;
				}
			}

		}//else
		
	}//while

	//if there is sequence of length less than Kmer in the end of chromosome, handle it as unique sequence
	if ((repeat2unique_flag == 1) && (pos >= seq.size()-Kmer+1) && (pos < seq.size()))
	{
		useq_count++;
		useq_len += seq.size() - pos;
		unique_outfile << ">" << useq_count << "\t" << pos << "\t" << seq.size()-1
			<< "\t" << seq.size() - pos << "\t" << id << "\n" << seq.substr(pos) << "\n";
	}

}

int main(int argc, char* argv[])
{
	if (argc < 2) { Usage(); return 0; }

	//handle command line
	int c;
	while ((c=getopt(argc, argv, "k:m:l:h")) != -1)
	{
		switch(c)
		{
			case 'k':	Kmer = atoi(optarg); break;
			case 'm':	init_size = atoi(optarg); break;
			case 'l':	load_factor = atof(optarg); break;
			case 'h':	Usage(); return 0;
			default:	Usage(); return 0;
		}
	}

	if (Kmer < 1 || Kmer > 60)
	{
		cerr << "Kmer=" << Kmer << " is out of range: [1, 60]\n";
		exit(1);
	}

	string genome_file = argv[optind++];
	ifstream infile(genome_file.c_str());
	if (!infile) { cerr << "Error: can't open input file: " << genome_file << endl;	exit(1); }


	//set log output file
	string log_output = genome_file + ".log";
	log_outfile.open(log_output.c_str());
	if (!log_outfile) { cerr << "Error: can't open log file: " << log_output << endl; exit(1); }

	//set kmer output file
	string kmer_output = genome_file + ".kmerFreqStat";
	kmer_outfile.open(kmer_output.c_str());
	if (!kmer_outfile) { log_outfile << "Error: can't open kmer file: " << kmer_output << endl; exit(1); }

	string kmerFreq_output = genome_file + ".kmerFreq";
	kmerFreq_outfile.open(kmerFreq_output.c_str());
	if (!kmerFreq_outfile) { log_outfile << "Error: can't open kmerSeq file: " << kmerFreq_output << endl; exit(1); }

	//set unique sequence output file
	string unique_output = genome_file + ".unique";
	unique_outfile.open(unique_output.c_str());
	if (!unique_outfile) { log_outfile << "Error: can't open output file: " << unique_output << endl; exit(1); }

	//set repeat sequence output file
	string repeat_output = genome_file + ".repeat";
	repeat_outfile.open(repeat_output.c_str());
	if (!repeat_outfile) { log_outfile << "Error: can't open output file: " << repeat_output << endl; exit(1); }

	log_outfile << argv[0] <<" -k " << Kmer << " -m " << init_size << " -l " << load_factor << " " << genome_file << "\n\n";

	start_time = time(NULL);

	//initialize hash table
	HashSet* Hash = init_hashset(init_size, load_factor);

	string line;

	//check the format of input file
	getline(infile, line, '\n');
	if (line[0] != '>') { log_outfile << "Sequence should be in *.fa format.\n"; exit(1); }

	//count number of chromosomes
	int chr_sum = 1;
	while (getline(infile, line, '>'))
	{
		chr_sum++;
	}
	infile.clear();
	chr_sum--;

	log_outfile << "There are " << chr_sum << " chromosomes.\n";

	string *chr_seq = new string[chr_sum];
	string *chr_id = new string[chr_sum];

	if (chr_seq == NULL || chr_id == NULL)
	{
		cerr << "Not enough memeory!\n";
		exit(3);
	}

	int current_chr = -1;
	string::size_type total_length = 0;

	//read input file
	infile.seekg(0, ios_base::beg); //re-locate the pointer to file's start position
	while (getline(infile, line, '\n'))
	{
		if (line[0] == '>')
		{
			if ((current_chr >= 0) && (chr_seq[current_chr].size() >= Kmer))
			{
				total_length += chr_seq[current_chr].size();

				//fill hash table
				fill_hash_table(Hash, chr_seq[current_chr]);
			}

			current_chr++;
			chr_id[current_chr] = line.substr(1);

			continue;
		}

		chr_seq[current_chr] += line;
	}
	infile.close();

	//handle last chromosome
	if ((current_chr >= 0) && (chr_seq[current_chr].size() >= Kmer))
	{
		total_length += chr_seq[current_chr].size();

		//fill hash table
		fill_hash_table(Hash, chr_seq[current_chr]);
	}

	log_outfile << "Total length is " << total_length << ".\n";

	log_outfile << "Count Kmer frequency done!\nThere are " << Hash->count << " different Kmers." << endl;


	//print depth of each appeared frequency, frequency >= 255 will be treated as 255
	print_hashset(Hash, kmer_outfile, kmerFreq_outfile);
	log_outfile << "Print kmer table done!" << endl;

	//get unique and repeat sequences and output them
	for (int i=0; i<chr_sum; i++)
	{
		if (chr_seq[i].size() >= Kmer)
		{
			fish_sequence(Hash, chr_seq[i], chr_id[i]);
		}
	}
	log_outfile << "Get uinique and repeat sequence done!" << endl;

	// free memory
	free_hash(Hash);
	delete []chr_seq;
	delete []chr_id;

	log_outfile << "There are " << useq_count << " unique sequences, " << rseq_count << " repeat sequences.\n"
			<< "The maximum length of unique sequence is " << max_unique_len 
			<< " , maximum length of repeat sequence is " << max_repeat_len << ".\n" 
			<< "Total length of unique sequences is " << useq_len << ".\n"
			<< "Total length of repeat sequences is " << rseq_len << ".\n" 
			<< "Total used time: " << time(NULL) - start_time << ".\nAll done!\n";
}

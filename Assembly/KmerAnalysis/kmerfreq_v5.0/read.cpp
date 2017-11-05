#include "read.h"
#include "threads.h"

extern int bufferNum;
extern string** readBuffer;
extern uint64_t* readNum;
extern uint64_t** indexArr;
extern uint64_t* kmerNum;
extern uint64_t maxReadNum;
extern int maxReadLen;

int check_seq (string &seq)
{	
	int is_good = 1;
	for (int i = 0; i < seq.size(); i++)
	{	if (alphabet[seq[i]] == 4)
		{	is_good = 0;
			break;
		}
	}
	return is_good;
}

bool ReadFilter(string &read_seq)
{
	if(Read_len > 0 && read_seq.size() < Read_len)
		return 0;
	return check_seq(read_seq);
}

void Loadfile(string &file, bit8_t* &freq)
{
        ifstream infile ( file.c_str() );

        if ( ! infile )
        {       cerr << "fail to open input file" << file << endl;
        }

        cerr<<"the operating file is: "<<file<<endl;
        //parse one reads at a time

        string seq;
        getline( infile, seq, '\n' );

        int type=0;
        if (seq[0] == '@'){type=1;} //the operating file is fq type.

	int buffId = 0;
        while ( getline( infile, seq, '\n' ) ) //read sequence
        {
                total_read++;

                if (! ReadFilter(seq))
		{
			if (type == 1)
			{
				getline(infile, seq, '\n');
				getline(infile, seq, '\n');
				getline(infile, seq, '\n');
			}
			else
			{
				getline(infile, seq, '\n');
			}
                        continue;
		}

                int real_len=0;
                if(Read_len == 1) {
			if (seq.size() < Delete)
			{
				if (type == 1)
				{
					getline(infile, seq, '\n');
					getline(infile, seq, '\n');
					getline(infile, seq, '\n');
				}
				else
				{
					getline(infile, seq, '\n');
				}
				continue;
			}
                        real_len=seq.size()-Delete;
                }else{
                        real_len=Read_len;
                }

		real_len = (real_len<maxReadLen) ? real_len : maxReadLen;
		if (real_len < Ahead+Kmer)
		{
			if (type == 1)
			{
				getline(infile, seq, '\n');
				getline(infile, seq, '\n');
				getline(infile, seq, '\n');
			}
			else
			{
				getline(infile, seq, '\n');
			}
			continue;
		}

		num_reads++;

                used_base+=real_len-Ahead;

		readBuffer[buffId][readNum[buffId]] = seq.substr(Ahead, real_len-Ahead);
		indexArr[buffId][readNum[buffId]] = kmerNum[buffId];
		kmerNum[buffId] += real_len-Ahead-Kmer+1;
		readNum[buffId]++;

		if (readNum[buffId] == maxReadNum)
		{
			num_kmers += kmerNum[buffId];

			SendSignal(buffId, 1);

			buffId = GetEmptyBuff();

			readNum[buffId] = 0;
			kmerNum[buffId] = 0;

		}

                if (num_reads % 10000000 == 0)
                {
                        cerr << "processed reads " << num_reads << "\n";
                }
                if((genome > 0) && (used_base >= genome)) //if set the total base
                {
                        cerr<<"Used total base: "<<genome<<endl;
                        break;
                }

                if(type == 1) //fq files
                {
                        getline( infile, seq, '\n' );
                        getline( infile, seq, '\n' );
                        getline( infile, seq, '\n' );
                }else{
                        getline( infile, seq, '\n' );
                }
        }
        infile.close();


	if (readNum[buffId] > 0)
	{
		num_kmers += kmerNum[buffId];

		SendSignal(buffId, 1);
	}

	WaitForThreads();

	for (int i=0; i<bufferNum; i++)
	{
		readNum[i] = 0;
		kmerNum[i] = 0;
	}
}

void getOtherLine_GZ(FILE *fp, int type)
{
	char buff[1024];

	fgets(buff, 1024, fp);

	if (type == 0)
	{
		return;
	}

	fgets(buff, 1024, fp);
	fgets(buff, 1024, fp);
}

void LoadGzfile(string &file, bit8_t* &freq)
{
	char *cmd = new char[file.size()+20];
	sprintf(cmd, "gzip -dc %s", file.c_str());

	FILE *fp = popen(cmd, "r");

/*
        gzFile zip;
        zip=gzopen (file.c_str(), "r");
        int c;
*/
        string seq;
        int line_num=0, type = -1;
        cerr<<"read file: "<<file<<endl;

	char buff[1024];
	char *line=NULL;

	fgets(buff, 1024, fp);
	seq = buff;

	if (seq[0] == '@')
	{
		type = 1;
	}
	else if (seq[0] == '>')
	{
		type = 0;
	}
	else
	{
		cerr << "Read file should be fastq or fasta format!\n";
		exit(1);
	}

	int buffId = 0;

/*
	line = gzgets(zip, buff, 1024);
	seq = line;
	if (line[0] == '@')
	{
		type = 1;
	}
*/
//	while (!gzeof(zip))
	while (fgets(buff, 1024, fp))
	{
//		line = gzgets(zip, buff, 1024);
//		seq = line;
		seq = buff;
		seq = seq.substr(0, seq.size()-1);

		total_read++;

                if (! ReadFilter(seq))
		{
			getOtherLine_GZ(fp, type);
/*
			if (type == 1)
			{

				gzgets(zip, buff, 1024);
				gzgets(zip, buff, 1024);
				gzgets(zip, buff, 1024);
			}
			else
			{
				gzgets(zip, buff, 1024);
			}
*/
                        continue;
		}

                int real_len=0;
                if(Read_len == 1) {
			if (seq.size() < Delete)
			{
				getOtherLine_GZ(fp, type);
/*
				if (type == 1)
				{
					gzgets(zip, buff, 1024);
					gzgets(zip, buff, 1024);
					gzgets(zip, buff, 1024);
				}
				else
				{
					gzgets(zip, buff, 1024);
				}
*/
				continue;
			}
                        real_len=seq.size()-Delete;
                }else{
                        real_len=Read_len;
                }

		real_len = (real_len<maxReadLen) ? real_len : maxReadLen;
		if (real_len < Ahead+Kmer)
		{
			getOtherLine_GZ(fp, type);
/*
			if (type == 1)
			{
				gzgets(zip, buff, 1024);
				gzgets(zip, buff, 1024);
				gzgets(zip, buff, 1024);
			}
			else
			{
				gzgets(zip, buff, 1024);
			}
*/
			continue;
		}

		num_reads++;

                used_base+=real_len-Ahead;

		readBuffer[buffId][readNum[buffId]] = seq.substr(Ahead, real_len-Ahead);
		indexArr[buffId][readNum[buffId]] = kmerNum[buffId];
		kmerNum[buffId] += real_len-Ahead-Kmer+1;
		readNum[buffId]++;

		if (readNum[buffId] == maxReadNum)
		{
			num_kmers += kmerNum[buffId];

			SendSignal(buffId, 1);

			buffId = GetEmptyBuff();

			readNum[buffId] = 0;
			kmerNum[buffId] = 0;

		}

                if (num_reads % 10000000 == 0)
                {
                        cerr << "processed reads " << num_reads << "\n";
                }
                if((genome > 0) && (used_base >= genome)) //if set the total base
                {
                        cerr<<"Used total base: "<<genome<<endl;
                        break;
                }

		getOtherLine_GZ(fp, type);
/*
                if(type == 1) //fq files
                {
			gzgets(zip, buff, 1024);
			gzgets(zip, buff, 1024);
			gzgets(zip, buff, 1024);
                }else{
			gzgets(zip, buff, 1024);
                }
*/		
	}

	if (readNum[buffId] > 0)
	{
		num_kmers += kmerNum[buffId];

		SendSignal(buffId, 1);
	}
//       gzclose(zip);

	pclose(fp);
	delete cmd;

	WaitForThreads();

	for (int i=0; i<bufferNum; i++)
	{
		readNum[i] = 0;
		kmerNum[i] = 0;
	}
}

void LoadList(string &files, bit8_t* &freq)
{
        ifstream list (files.c_str());
        if(!list)
        {
                cerr<<"fail to open input file: "<<files<<endl;
        }
        string file;
        while(getline (list, file,'\n'))
        {
                cerr<<"File: "<<file<<endl;

		string suffix = file.substr(file.size()-3);
		if (suffix == ".gz")
//		if(file.find(".gz") != string::npos)
                {
                        LoadGzfile(file, freq);
                }else{
                        Loadfile(file, freq);
                }
        }
        list.close();
        cerr<<"Load All files done!\n";
}

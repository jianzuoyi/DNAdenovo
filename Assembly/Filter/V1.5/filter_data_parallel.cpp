#include <iostream>
#include <string>
#include <cstdlib>
#include <stdint.h>
#include <pthread.h>

#include "filter_low_quality.h"
#include "filter_adapter.h"
#include "filter_small_size.h"
#include "gzstream.h"

using namespace std;

int Q_SHIFT = 64;

int Start_trim1 = 0;
int End_trim1 = 0;
int Start_trim2 = 0;
int End_trim2 = 0;

float N_rate = 0.1;
int N_num = -1;
int Qual_rate=40;

int Do_filter_Ns = 1;
int Do_low_qual = 1;
int Do_adapter = 0;
int Do_small_size = 0;
int insert_size=500;
int buffNum=2000000;

string* read1ID = NULL;
string* read1Seq = NULL;	//store reads
string* read1Quality = NULL;
string* read2ID = NULL;
string* read2Seq = NULL;
string* read2Quality = NULL;
int* IsDelete = NULL;
int* gotReads = NULL;
int threadNum = 8;

igzstream infile1,infile2;
ofstream stat,clean1,clean2;

int Read_len1 = 0,Read_len2 = 0,Usable_len1 = 0,Usable_len2 = 0;;
uint64_t Raw_reads = 0,Raw_bases = 0,Read1_GC = 0,Read2_GC = 0,Read1_Q20 = 0,Read2_Q20 = 0;
uint64_t Low_qual_num = 0, Ns_num = 0, Adapter_num = 0, Small_num = 0;
uint64_t readNum;

void* thread_filter_data(void* threadId_p);
void usage ();

int main(int argc, char *argv[])
{
	if (argc<5)
	{
		usage();
	}
	int c;
	while ((c=getopt(argc,argv,"a:b:c:d:q:m:B:l:w:t:yzh")) != -1)
	{
		switch (c)
		{
			case 'a' : Start_trim1 = atoi(optarg); break;
			case 'b' : End_trim1 = atoi(optarg); break;
			case 'c' : Start_trim2 = atoi(optarg); break;
			case 'd' : End_trim2 = atoi(optarg); break;
			case 'B' : Qual_rate = atoi(optarg);break;
			case 'l' : insert_size = atoi(optarg);break;
			case 'q' : Q_SHIFT = atoi(optarg);break;
			case 'm' : buffNum = atoi(optarg);break;
			case 't' : threadNum = atoi(optarg);break;
			case 'w' : N_num = atoi(optarg); break;
			case 'y' : Do_adapter = 1; break;
			case 'z' : Do_small_size =1; break;
			case 'h' : usage(); break;
			default  : cout<<"error:"<<(char)c<<endl;usage();
		
		}
	}
	
	if(N_num>=0){N_rate=float(N_num)/100;}

	string seq1_file = argv[optind++]; //optind, argv[optind++]顺序指向非option的参数
	string seq2_file = argv[optind++];
	string stat_file = argv[optind++];
	string seq1_clean = argv[optind++];
	string seq2_clean = argv[optind++];

	infile1.open(seq1_file.c_str());
	if(!infile1)
	{
		cerr<<"fail to open input file : "<<seq1_file<<endl;
		exit(1);
	}
	
	infile2.open(seq2_file.c_str());
	if(!infile2)
	{
		cerr<<"fail to open input file : "<<seq2_file<<endl;
		exit(1);
	}
	
	stat.open(stat_file.c_str());
	if(!stat)
	{
		cerr<<"fail to create stat file : "<<stat_file<<endl;
		exit(1);
	}
	
	clean1.open(seq1_clean.c_str());
	if(!clean1)
	{
		cerr<<"fail to create read1 clean file : "<<seq1_clean<<endl;
		exit(1);
	}
	
	clean2.open(seq2_clean.c_str());
	if(!clean2)
	{
		cerr<<"fail to create read2 clean file : "<<seq2_clean<<endl;
		exit(1);
	}
	
	read1ID = new string[buffNum];
	read1Seq = new string[buffNum];
	read1Quality = new string[buffNum];
	read2ID = new string[buffNum];
	read2Seq = new string[buffNum];
	read2Quality = new string[buffNum];
	IsDelete = new int[buffNum];
	gotReads = new int[threadNum];
	string qid;
	
	for(int i=0;i<buffNum;i++)
		IsDelete[i] = 0;

	pthread_t *pthread = new pthread_t[threadNum];
	int *pthreadId = new int[threadNum];
	
	for (int i=0; i<threadNum; i++)
	{	
		gotReads[i] = 0;
		pthreadId[i] = i;
		pthread_create((pthread+i), NULL, thread_filter_data, (void*)(pthreadId+i));
	}
	cerr << threadNum << " threads creation done!\n" << endl;
	
	
// 获得截取后的有效reads长度以及判断是否进行small_insert
	
	string textLine,seq1="",seq2="";
	igzstream in1(seq1_file.c_str());
	if(!in1)
	{
		cerr<<"fail to open input file : "<<seq1_file<<endl;
		exit(1);
	}
	for(int i=0;i<100 && getline(in1,textLine,'\n');i++)
	{
		if(textLine[0] == '@')
		{
			getline(in1,seq1,'\n');
			break;
		}
	}
	igzstream in2(seq2_file.c_str());
	if(!in2)
	{
		cerr<<"fail to open input file : "<<seq2_file<<endl;
		exit(1);
	}
	for(int i=0;i<100 && getline(in2,textLine,'\n');i++)
	{
		if(textLine[0]=='@')
		{
			getline(in2,seq2,'\n');
			break;
		}
	}
	if(seq1.size()==0 || seq2.size()==0)
	{
		cerr<<"please make sure the input files are fastq format !!!"<<endl;
		exit(2);
	}
	
	Read_len1 = seq1.length();
	Usable_len1 = Read_len1 - Start_trim1 - End_trim1;
	Read_len2 = seq2.length();
	Usable_len2 = Read_len2 - Start_trim2 - End_trim2;
	if(Read_len1+Read_len2+30 >insert_size)
		Do_small_size=0;
	in1.close();
	in2.close();
	
	while(1)
	{
		readNum = 0;
		while(readNum<buffNum && getline(infile1,read1ID[readNum],'\n') && getline(infile2,read2ID[readNum],'\n'))
		{
//cout<<readNum<<'\t'<<read1ID[readNum]<<endl;
			getline(infile1,read1Seq[readNum],'\n');
			getline(infile1,qid,'\n');
			getline(infile1,read1Quality[readNum],'\n');
			
			//stat GC content and Q20 of read1
			for(int i=0;i<Read_len1;i++)
			{
				switch(read1Seq[readNum][i])
				{
					case 'G':
					case 'C': Read1_GC++;break;
				}
				if(read1Quality[readNum][i] - Q_SHIFT >= 20)
					Read1_Q20++;
			}
			
			getline(infile2,read2Seq[readNum],'\n');
			getline(infile2,qid,'\n');
			getline(infile2,read2Quality[readNum],'\n');
			
			Raw_reads++;
			Raw_bases += Read_len1;
			Raw_bases += Read_len2;
			
			//stat GC content and Q20 of read1
			for(int i=0;i<Read_len2;i++)
			{
				switch(read2Seq[readNum][i])
				{
					case 'G':
					case 'C': Read2_GC++;break;
				}
				if(read2Quality[readNum][i]-Q_SHIFT >= 20)
					Read2_Q20++;
			}
			readNum++;
		}
		
		//线程运行，直到Loop_id等于readNum, 然后进入子线程等待状态
		//放行主线程，进入下一个循环，读数据
		for (int i=0; i<threadNum; i++)
		{	gotReads[i] = 1;
		}
		
		//等子线程，直到gotReads都是0，说明任务完成，可以继续读下一批数据
		while (1)
		{
			usleep(1);
			int i=0;
			for (; i<threadNum; i++)
			{	if (gotReads[i] == 1)
				{	break;
				}
			}
			if (i == threadNum)
			{	break;
			}
		}
		
		//输出结果
		for(int i=0;i<readNum;i++)
		{
			if(IsDelete[i] == 0)
			{
				clean1<<read1ID[i]<<'\n'<<read1Seq[i]<<"\n+\n"<<read1Quality[i]<<endl;
				clean2<<read2ID[i]<<'\n'<<read2Seq[i]<<"\n+\n"<<read2Quality[i]<<endl;
			}
			else
				IsDelete[i] = 0;
		}

		//如果已经到了文件尾，则终止全部子线程，并且退出读文件的循环
		if (readNum < buffNum)
		{	for (int i=0; i<threadNum; i++)
			{	gotReads[i] = 2;
			}
			break; //读到文件尾退出
		}
	}
	
	//等待全部子线程结束
	for (int i=0; i<threadNum; i++)
	{
		pthread_join(pthread[i], NULL);
	}

	infile1.close();
	infile2.close();
	clean1.close();
	clean2.close();
	
	//output stat. information
	stat.precision(3);
	double Q201 = 2*100* ((double)Read1_Q20)/Raw_bases;
	double Q202 = 2*100* ((double)Read2_Q20)/Raw_bases;
	double GC1 = 2*100* ((double)Read1_GC)/Raw_bases;
	double GC2 = 2*100* ((double)Read2_GC)/Raw_bases;
	stat << "Raw_reads\tRaw_len\tRaw_bases\tRead1_Q20\tRead2_Q20\tRead1_GC\tRead2_GC\t"
			<<"Ns_num\tLow_qual\tAdapter\tSmall\tUsable_len1\tUsable_len2" << endl;
	stat <<Raw_reads<<"\t"<<Read_len1<<"_"<<Read_len2<<"\t"<<Raw_bases<<"\t"<<Q201<<"\t"<<Q202<<"\t"<<GC1<<"\t"<<GC2<<"\t" 
			<<Ns_num<<"\t"<<Low_qual_num<<"\t"<<Adapter_num<<"\t"<<Small_num<<"\t"<<Usable_len1<<"\t"
			<<Usable_len2<<endl;
	stat.close();
	
}

void usage ()
{
	cout << "\nUsage: filter_data_parallel [options] <read_1.fq> <read_2.fq> <raw_reads_stat> <read_1.clean> <read_2.clean>\n"
			<< "  -a <int>  trimed length at 5' end  of read1, default " << Start_trim1 << "\n"
			<< "  -b <int>  trimed length at 3' end of read1, default " << End_trim1 << "\n"
			<< "  -c <int>  trimed length at 5' end of read2, default " << Start_trim2 << "\n"
			<< "  -d <int>  trimed length at 3' end of read2, default " <<End_trim2 << "\n"
			<< "  -q <int>  the quality shift value 64 or 33, default " <<Q_SHIFT <<"\n"
			<< "  -m <int>  the reads pair number in buffer,default " <<buffNum <<"\n"
			<< "  -t <int>  run the program in multiple thread mode, default=" << threadNum << endl
			<< "  -B <int>  filter reads with many low quality bases,set a cutoff, default " << Qual_rate << "\n"
			<< "  -l <int>  library insert size, default " <<insert_size << "\n"
			<< "  -w <int>  filter reads with >X percent base is N,set a cutoff, default 10" << endl
			<< "  -y        filter reads with adapter " << endl
			<< "  -z        filter reads with small size " << endl
			<< "  -h        output help information\n" << endl;
	exit(1);
}

void* thread_filter_data(void* threadId_p)
{
	int threadId = *((int*)threadId_p);
	while (1)
	{
		usleep(1);

		if (gotReads[threadId] == 1)
		{
			for (int i=0; i<readNum; i++)
			{
				if (i%threadNum == threadId)
				{
					if ( Do_adapter && filter_adapter(read1ID[i],read1Seq[i],read2ID[i],read2Seq[i]) )//filter adapter
					{
						__sync_add_and_fetch(&Adapter_num,1);
						IsDelete[i]=1;
						continue;
					}
					if ( Start_trim1 > 0 || End_trim1 >0 || Start_trim2 >0 ||End_trim2>0)
						trim(read1Seq[i], read2Seq[i], read1Quality[i], read2Quality[i], Start_trim1, End_trim1, Start_trim2, End_trim2);//trim reads
					if ( Do_filter_Ns && filter_Ns(read1Seq[i],read2Seq[i],N_rate))//filter Ns and polyA
					{
						__sync_add_and_fetch(&Ns_num,1);
						IsDelete[i]=1;
					}
					else if ( Do_low_qual && filter_low_qual(read1Quality[i],read2Quality[i],Qual_rate) ) //filter_low_qual
					{
						__sync_add_and_fetch(&Low_qual_num,1);
						IsDelete[i]=1;
					}
					else if ( Do_small_size && filter_small_size(read1Seq[i],read2Seq[i]) )//filter small insert size
					{
						__sync_add_and_fetch(&Small_num,1);
						IsDelete[i]=1;
					}
				}
			}

			gotReads[threadId] = 0;
		}
		else if (gotReads[threadId] == 2)
		{
			return NULL;
		}
	}
}

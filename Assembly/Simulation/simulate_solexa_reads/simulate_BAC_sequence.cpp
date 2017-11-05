#include <iostream>
#include <string>
#include <fstream>
#include <vector>
using namespace std;

string input;
char *output;
int BAC_length=100000;
int BAC_number=100;
int BAC_group=1;

const char *MAKE_TIME="2010-05-04";
const char *VERSION="1.0";
const char *AUTHOR="lujianliang";
const char *CONTACT="lujianliang@genomics.org.cn";

void Usage(){
	cout<<"Description:"<<endl;
	cout<<endl<<"\tIt's a program simulate BAC sequence or fosmid sequence."<<endl;
	cout<<"\tThe length of simulated sequence are all the same and the length can be set by user,"<<endl;
	cout<<"\tif the reference sequence is multiple,each sequence generate BAC or fosmid sequence"<<endl;
	cout<<"\tin proportional with their length.the parameter -g is used to simulate sequence of"<<endl;
	cout<<"\tthe different index.the meanning of output file name is ,for example,BAC_100000_100_1_.fa,"<<endl;
	cout<<"\tthe 'BAC' is output prefix,the '100000' is the length of BAC sequence,the '100' is the BAC"<<endl;
	cout<<"\tnumber to be simulate,and the last '1' is the group number that set by -g.The meanning of"<<endl;
	cout<<"\tsequence id is,for example,I_71610_100000_1,the 'I' is the reference sequence id,the '71610'"<<endl;
	cout<<"\tis the start position of simulated sequence,the start position of the reference sequence is 1."<<endl;
	cout<<"\tthe '100000' is the length of simulated sequence,and the last '1' is the sequence id."<<endl;
	cout<<"\tCompile Data:\t"<<MAKE_TIME<<endl;
	cout<<"\tAuthor:\t\t"<<AUTHOR<<endl;
	cout<<"\tVersion:\t"<<VERSION<<endl;
	cout<<"\tContact:\t"<<CONTACT<<endl;
	cout<<endl<<"Usage:\tsimulate_BAC_sequence [options]"<<endl;
	cout<<"\t-i 	<string>	input reference genome sequence *.fa"<<endl;
	cout<<"\t-l	<int>		input the length of BAC sequecne(default:100000)"<<endl;
	cout<<"\t-n	<int>		input the number of BAC sequence(default:100)"<<endl;
	cout<<"\t-g	<int>		input the group(index) number(default:1)"<<endl;
	cout<<"\t-o	<string>	output file prefix(default:BAC)"<<endl;
	cout<<"\t-h			output help information"<<endl;
	cout<<endl<<"Example:"<<endl;
	exit(-1);
}

void Getopt(int argc,char *argv[]){
	int c;
	while ((c=getopt(argc,argv,"i:l:n:g:o:h"))!=-1)
	{
		switch(c){
			case 'i': input=optarg;break;
			case 'l': BAC_length=atoi(optarg);break;
			case 'n': BAC_number=atoi(optarg);break;
			case 'g': BAC_group=atoi(optarg);break;
			case 'o': output=optarg;break;
			case 'h': Usage();break;
			default: Usage();
		}
	}
}

void get_ref(ifstream &inputf,vector<string> &sequence1,vector<string> &id1,vector<int> &chr_length1,long long &genome_length1){
	string line,id,seq;
	int len=0;
	while (getline(inputf,line,'\n'))
	{
		if (line[0]=='>')
		{
			if (seq!="")
			{	
				sequence1.push_back(seq);
				id1.push_back(id);
				len=seq.size();
				chr_length1.push_back(len);
				genome_length1+=len;
				seq="";
			}
			line.erase(0,1);
			int pos=line.find(" ");
			line=line.substr(0,pos);
			id=line;
		}else{
			seq+=line;
		}		
	}
	sequence1.push_back(seq);
	id1.push_back(id);
	len=seq.size();
	chr_length1.push_back(len);
	genome_length1+=len;
}

void get_BAC(vector<string> &sequence2,vector<string> &id2,vector<int> &chr_length2,long long genome_length2,char *output1){
	int num=0;int a;
	for (a=0;a<chr_length2.size();a++)
	{
		chr_length2[a]=(int)(BAC_number*(double)chr_length2[a]/(double)genome_length2);
		num+=chr_length2[a];
	}
	if (num<BAC_number)
	{
		chr_length2[--a]+=BAC_number-num;
	}
	num=0;
	while (BAC_group>0)
	{
		int num2=0;
		num++;
		char output2[500];
		sprintf(output2,"%s%s%d%s",output1,"_",num,".fa");
		ofstream out;
		out.open(output2);
		if (!out)
		{
			cout<<"Error: can't open file "<<output<<endl;
			exit(-1);
		}
		for (int b=0;b<sequence2.size();b++)
		{
			int size=sequence2[b].size();
			int seq_num=chr_length2[b];
			while (seq_num>0)
			{
				int pos=int (rand()%size);
				if (pos+BAC_length>size) continue;
				string BAC=sequence2[b].substr(pos,BAC_length);
				num2++;
				out<<">"<<id2[b]<<"_"<<pos+1<<"_"<<BAC_length<<"_"<<num2<<endl;
				out<<BAC<<endl;
				seq_num--;
			}
		}
		out.close();
		BAC_group--;
	}
}

int main(int argc, char *argv[])
{
	if (argc==1) Usage();
	Getopt(argc,argv);
	ifstream inputfile;
	inputfile.open(input.c_str());
	if (!inputfile)
	{
		cout<<"Error: can't open input file "<<input<<endl;
		exit(-1);
	}
	char output1[500];
	if (!output)
	{
		sprintf(output1,"%s%d%s%d","BAC_",BAC_length,"_",BAC_number);
	}else{
		sprintf(output1,"%s%s%d%s%d",output,"_",BAC_length,"_",BAC_number);
	}

	vector<string> sequence;
	vector<string> id;
	vector<int> chr_length;
	long long genome_length=0;
	get_ref(inputfile,sequence,id,chr_length,genome_length);
	get_BAC(sequence,id,chr_length,genome_length,output1);
	return 0;
}



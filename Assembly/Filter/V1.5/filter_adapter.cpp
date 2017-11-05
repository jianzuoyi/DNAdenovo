#include "filter_adapter.h"

using namespace std;

int align(string seq,string adpt,int& seq_start,int& seq_end,int& adpt_start,int& adpt_end);
string reverse(string seq);

int filter_adapter( string& id1, string& seq1, string& id2, string& seq2)
{
	string adpt1="AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG";//pDNA-3+
	string adpt2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT";//pDNA-5-
	int seq1_start=0,seq1_end=0,adpt1_start=0,adpt1_end=0;
	int seq2_start=0,seq2_end=0,adpt2_start=0,adpt2_end=0;
	int min_score=38;
	int score1 = align(seq1,adpt1,seq1_start,seq1_end,adpt1_start,adpt1_end);
	int score2 = align(seq2,adpt2,seq2_start,seq2_end,adpt2_start,adpt2_end);
	if(score1>=min_score || score2>=min_score)
	{
		int sub_len = seq1_start-adpt1_start > seq2_start-adpt2_start ? seq2_start-adpt2_start : seq1_start-adpt1_start;
		if(sub_len >= 10)
		{
			string seq22 = reverse(seq2);
			int mis_num = 0;
			int max_mis_num = int(sub_len * 0.1);
			for(int i=0,j=seq22.length()-sub_len;i<sub_len;i++,j++)
			{
				mis_num += !(!(seq1[i]-seq22[j]));
				if(mis_num > max_mis_num)
					return 0;
			}
//			cout<<id1<<'\t'<<sub_len<<'\t'<<score1<<'\t'<<seq1_start<<'\t'<<seq1_end<<'\t'<<adpt1_start<<'\t'<<adpt1_end<<endl;
//			cout<<id2<<'\t'<<sub_len<<'\t'<<score2<<'\t'<<seq2_start<<'\t'<<seq2_end<<'\t'<<adpt2_start<<'\t'<<adpt2_end<<endl;
			return 1;
		}
		else if(sub_len>=0)
		{
			int adpt_len = adpt1.length();
			int mis_num1=0,mis_num2=0,max_mis_num = int(adpt_len * 0.1);
			for(int i=sub_len,j=0;j<adpt_len;i++,j++)
			{
				mis_num1 += !(!(seq1[i]-adpt1[j]));
				mis_num2 += !(!(seq2[i]-adpt2[j]));
				if(mis_num1>max_mis_num && mis_num2>max_mis_num)
					return 0;
			}
//			cout<<id1<<'\t'<<sub_len<<'\t'<<score1<<'\t'<<seq1_start<<'\t'<<seq1_end<<'\t'<<adpt1_start<<'\t'<<adpt1_end<<endl;
//			cout<<id2<<'\t'<<sub_len<<'\t'<<score2<<'\t'<<seq2_start<<'\t'<<seq2_end<<'\t'<<adpt2_start<<'\t'<<adpt2_end<<endl;
			return 1;
		}
		else
		{
			int adpt_len = adpt1.length();
			int mis_num1=0,mis_num2=0;
			int match1 = seq1_end+adpt_len-adpt1_end,match2=seq2_end+adpt_len-adpt2_end;
			if(match1<10 || match2<10)
				return 0;
			int max_mis_num1=int(0.1*(seq1_end+adpt_len-adpt1_end)),max_mis_num2=int(0.1*(seq2_end+adpt_len-adpt2_end));
			for(int i=match1-1,j=adpt_len-1;i>=0;i--,j--)
			{
				mis_num1 += !(!(seq1[i]-adpt1[j]));
				if(mis_num1>max_mis_num1)
					return 0;
			}
			for(int i=match2-1,j=adpt_len-1;i>=0;i--,j--)
			{
				mis_num2+= !(!(seq2[i]-adpt2[j]));
				if(mis_num2>max_mis_num2)
					return 0;
			}
//			cout<<id1<<'\t'<<sub_len<<'\t'<<score1<<'\t'<<seq1_start<<'\t'<<seq1_end<<'\t'<<adpt1_start<<'\t'<<adpt1_end<<endl;
//			cout<<id2<<'\t'<<sub_len<<'\t'<<score2<<'\t'<<seq2_start<<'\t'<<seq2_end<<'\t'<<adpt2_start<<'\t'<<adpt2_end<<endl;
			return 1;
		}
	}
	return 0;
}


int align(string seq,string adpt,int& seq_start,int& seq_end,int& adpt_start,int& adpt_end)
{
	int seq_len = seq.length();
	int adpt_len = adpt.length();
	int j_size = adpt_len + 1;
	int *DPscore = new int[(seq_len+1)*(adpt_len+1)];
	int score[2]={-7,5};
	for(int i=0;i<=seq_len;i++)
		DPscore[i*j_size + 0]=0;
	for(int j=0;j<=adpt_len;j++)
		DPscore[0*j_size + j]=0;
	for(int i=1;i<=seq_len;i++)
		for(int j=1;j<=adpt_len;j++)
		{
			DPscore[i*j_size + j] = DPscore[(i-1)*j_size + (j-1)] + score[!(seq[i-1]-adpt[j-1])];
			if(DPscore[i*j_size + j]<0)
				DPscore[i*j_size + j] = 0;
		}
	int max_score = 0;
	for(int i=0;i<=seq_len;i++)
		for(int j=0;j<=adpt_len;j++)
			if(DPscore[i*j_size + j]>max_score)
			{
				max_score = DPscore[i*j_size + j];
				seq_end = i;
				adpt_end = j;
			}
	for(int i=seq_end,j=adpt_end;DPscore[i*j_size+j]>0;i--,j--)
	{
		seq_start=i;
		adpt_start=j;
	}
	delete[] DPscore;
	return max_score;

}

string reverse(string seq)
{
	int length=seq.length();
	string str;
	for(int i=length-1;i>=0;i--)
	{
		if (seq[i]=='A')
			str.push_back('T');
		else if(seq[i]=='T')
			str.push_back('A');
		else if(seq[i]=='C')
			str.push_back('G');
		else if(seq[i]=='G')
			str.push_back('C');	
		else
			str.push_back('N');
	}
	return str;
}

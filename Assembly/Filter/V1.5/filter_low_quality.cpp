#include "filter_low_quality.h"

using namespace std;

void trim(string &seq1, string &seq2, string &qual1, string &qual2, int start_trim1, int end_trim1, int start_trim2, int end_trim2)
{
	int len1 = seq1.length();
	int len2 = seq2.length();
	if ( start_trim1 + end_trim1 >=len1 || start_trim2 + end_trim2 >=len2 )
	{
		cerr << " Trimed str longer than read length "<<endl;
	}
	string s1="", q1="", s2="", q2="";
	int i=0, j=0;
	if (start_trim1 >0 || end_trim1 >0){
		for ( i=start_trim1;i<len1-end_trim1;i++ )
		{
			s1+=seq1[i];
			q1+=qual1[i];
		}
		seq1 = s1;
		qual1 = q1;
	}else{
		//
	}
	if (start_trim2 > 0 || end_trim2 > 0){
		for ( j=start_trim2;j<len2-end_trim2;j++ )
		{
			s2+=seq2[j];
			q2+=qual2[j];
		}
		seq2 = s2;
		qual2 = q2;
	}
}

int filter_Ns(string &seq1, string &seq2, float N_rate)
{
	float cutoff1 = seq1.length() * N_rate;
	float cutoff2 = seq2.length() * N_rate;
	int n1=0, n2=0;
	int flag=0;
	for ( int i=0;i<seq1.length();i++ )
	{
		if ( seq1[i] == 'N' )
		{
			n1++;
			if ( n1 > cutoff1 )
			{
				return 1;
			}
		}else if ( seq1[i] != 'A' )
		{
			flag=1;
		}
	}
	for ( int j=0;j<seq2.length();j++ )
	{
		if ( seq2[j] == 'N' )
		{
			n2++;
			if ( n2 > cutoff2 )
			{
				return 1;
			}
		}else if ( seq2[j] != 'A' )
		{
			flag=1;
		}
	}
	if ( flag)
	{
		flag=0;
		return 0;
	}else{
		return 1;
	}
}

int filter_low_qual(string &qual1, string &qual2, int Qual_rate)
{
	int Q_shift=-64;
	int len1 = qual1.length();
	int len2 = qual2.length();
	int qlen1 = Qual_rate;
	int qlen2 = Qual_rate;
	int q1=0,q2=0;
	for ( int i=0;i<len1;i++ )
	{
		if ( qual1[i]+Q_shift<=7)
		{
			q1++;
			if ( q1 >= qlen1 )
			{
				return 1;
			}
		}
	}
	for (int j=0;j<len2;j++ )
	{
		if ( qual2[j]+Q_shift<=7 )
		{
			q2++;
			if ( q2 >= qlen2 )
			{
				return 1;
			}
		}
	}
	return 0;
}

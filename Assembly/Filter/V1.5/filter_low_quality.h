#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std;

int filter_Ns(string &seq1, string &seq2, float N_rate);
int filter_low_qual(string &qual1, string &qual2, int Qual_rate);
void trim(string &seq1, string &seq2, string &qual1, string &qual2, int start_trim1, int end_trim1,int start_trim2, int end_trim2);

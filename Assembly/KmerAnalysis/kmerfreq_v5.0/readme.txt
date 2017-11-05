To use the program:
1. prepare the pathway of fq files.
2. qsub the work.sh

if k=17, then use 16G.
----------------------------------------------------------
Result:

creat: *.freq; *.log.
*.freq:
depth	frequence	frequency	frequence*depth

*.log:
tail -4 *.log
show the statistic table.

you can use many software to plot.
here we provide the gnuplot.

3. copy the creat_gnu.pl to your dir.
4. perl creat_gnu.pl >kmer.gnu;
5. gnuplot kmer.gnu;

creat: *.kmer.freq.png  
6. display *.kmer.freq.png;

The raw picture may be ugly, so you need to learn more about gnuplot and change parametors in kmer.gnu. Or, you can use other plot software.
------------------------------------------------------------
Attention:
1. Do not set -r and -d at the same time.
2. Do not set -k as a even number(2*int).
-----------------------------------------------------------
History:
The first author of the program is Li Ruiqiang, and modified by Fan Wei from 2009-3-18 up to now.
After version kmerfreq 4.0, Li Zhenyu support many great functions.

-----------------------------------------------------------
Any questions, please contact: liubinghang@genomics.org.cn.

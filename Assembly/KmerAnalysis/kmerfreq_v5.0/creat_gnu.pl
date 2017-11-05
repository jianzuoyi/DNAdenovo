#!user/bin/perl

#used to creat kmerfreq.gnu to plot the kmer distribution using gnuplot.
$dir=` pwd `;
chomp($dir);
#$dir=shift;
@files=glob("$dir/*.freq");

#creat the basic of the gnu file
#print "set terminal postscript portrait color solid  \"Times-Roman\" 30;\n"; #for pdf 
print "set terminal png \n";#\"Times-Roman\" 30;\n"; #for pdf
print "set xrange [1:100];\n";#set yrange [0:5];\n";
print "set xlabel \"Depth\";set ylabel \"Frequency(%)\";\n";
print "set xtics 1,10,100;set ytics nomirror;set xtics nomirror;\n";
print "set size 2.0,1;\n";#set border lw 2;\n";
print "set key top right;set key box lw 2;\n";

foreach $file(@files){
	print "set output \"$file.png\";\n";
	print "plot \"$file\" using 1:3 title \"Kmer_frequency\" with lines lw 3;\n";
	print "plot \"$file\" using 1:4 title \"Product_frequency\" with lines lw 3;\n";	
}
#print "replot\n";
#print "set terminal x11;\n";

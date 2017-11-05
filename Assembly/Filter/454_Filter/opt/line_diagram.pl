#!/usr/bin/perl  -w

=head1 Name

 line_diagram.pl

=head1 Description

  to draw line diagram, also it can draw bar or dot chart

=head1 Version

 Author: Wenbin Liu, liuwenbin@genomics.org.cn
 Version: 1.1,  Date: 2010-11-29
 Version: 2.0,  Date: 2011-01-19
                Update: add option --barstroke,barfill,opacity,dotsel,nolsel,barsel,trate
 Version: 2.1,  Date: 2011-03-04
                Update: add option -h_title -size_h -border
 Version: 3.0,  Date: 2011-05-27
                Update: add option -fontfam -italic
 Version: 3.1,  Date: 2011-7-26
                Update: add option -windl,-sym_xy,'pxxx,pyyy', also .gz infile is allow

=head1 Usage

 perl line_diagram.pl  <infiles>  [--Option --help]  > out.svg
 infiles            files store Data for drawing, the data can store in more then one file
 1 about statistics:
 --fredb            to calculate frequence distribution data to draw figure
 --group <num>      group number for frequence calculated, when --fredb, default=50
 --windl <flo>      wind length for frequence stated, default stated acroding to --group
 --valgp            one value one group when --fredb, then --group failure
 --numberc          y-axis data not to be frequency but number when --fredb
 --ranky <str>      y-axis data source: file0:rank1,rank2..;filei:rankj,rankk.., default 0:2
 --rankx <str>      x-axis data source, form as --ranky, default 0:1, failure when --fredb
 --samex            all the y-axis stats data use the same x-axis data
 --row              the data store in rows not in ranks
 --signh            the head of rank not state date but the symbol
 --ignoz            ignore the data equal zero
 --plot <file>      output x,y plot date for darwing to specified file
 
 
 2 main draw option
 2.1 about paper size
 --width <num>      the length of the x axis,default=400 
 --height <num>     the length of the y axis,default=300
 --flank_x <num>    distance from x axis to edge of drawing paper, default=width/4
 --flank_y <num>    distance from y axis to edge of drawing paper, default=height/3
 
 2.2 about dot
 --dot              show the key points at the line of all the group
 --dotsel <str>     selected line show dots, '1' means the first line, failure when --dot
 --size_d <num>     the size of the dot, default = 0.008*width
 --onedot           use only sycle sign when -dot
 
 2.3 about line
 --noline           not show the line of all the group
 --nolsel <str>     selected line not show line, '1' means the first line, failure when --noline
 --linew <num>      the line-stroke width, default = 2
 --color <str>      the colors of the lines, splited by ',' or in a file, default auto
 --border           line can not outsize border
 
 2.4 about bar
 --bar              to draw the bar chart instead of line chart
 --barsel <str>     selected line draw the bar, '1' means the first line, failure when --bar
 --barfill <str>    the bar fill colors, default set auto by the process
 --barstroke <str>  the bar stroke colors, default the same as the fill color
 --opacity <str>    the percentage of bar stroke-opacity, default=100
 
 2.5 about symbol
 --symbol           show the symbol
 --signs <str>      the sign of each line, splited by ',' or in a file, default get from --signh
 --size_st <num>    the font-size of symbol text, default=0.035*width
 --size_sg <num>    the size of symbol signs dot, default = 0.008*width
 --syml <num>       symbol line length, default = 0.05*width
 --sym_xy           symbol start coordinate x,y in xoy-axis figure, default top left corner
                    it can also head with 'p', means percent of the figure, e.g: p0.5,p0.9.
 --sym_frame        add the symbol frame
 
 2.6 about title
 --h_title <str>    the head title, default no title
 --x_title <str>    the figure x-area title,default no title
 --y_title <str>    the figure y-area title,default no title
 --size_h <num>     the head title font-size, default = 0.04*width
 --size_xt <num>    the font-size of x-area title text, default=0.035*width
 --size_yt <num>    the font-size of y-area title text, default=0.08*height
 --fontfam <str>    font-family, default=Arial
 --italic           text with [] in title use font-sytle italic
 --golden           put the y-title center at the Golden-Point of y-axis, default middle
 
 2.7 about scale
 --rever            reverse the data turn according to x-area data
 --frame            show the frame
 --gridx            show x-axis vertical grid
 --gridy            show y-axis vertical grid
 --micrx            show x-axis micro scale
 --micry            show y-axis micro scale
 --micrf            show micro scale in frame
 --scalen           scale short line show entad
 --miceq <str>      {main,micrf scale line length}/scale-font-size, default '0.35,0.25';
 --size_xs <num>    the font-size of x-axis scale text,default=0.035*width
 --size_ys <num>    the font-size of y-area scale text,default=0.035*width
 --x_mun <str>      "min,unit,scale_number" show at x-axis, default base on input date
 --y_mun <str>      "min,unit,scale_number" show at y-axis, default base on input date
 --x_scale <str>    the real x-axis scale showed, splited by ',' or set infile, default auto
 --scalmid          show x-area scale at the middle position
 --scalx <num>      the time to zoom in the x-axis data, deafult not zoom
 --scaly <num>      the time to zoom in the y-axis data, deafult not zoom
 --prex <num>       x-axis scale precision option, default=4
 --prey <num>       y-axis scale precision option, default=3
 --trate <num>      the rate of max-yvalue to max-yscale, default=95
 --edgex            x-axis scale goto the edge of the max x_value
 --edgey            x-axis scale goto the edge of the max y_value
 --size_zo <flo>    to zoom in all the text forn-size, default=1
 
 
 3 about vice
 --vice             draw vice y-axis line, follow options can use with meaning as upper
                    --fredb2 --valgp2 --ranky2 --rankx2 --y_title2 --bar2 --y_mun2 --prey2
                    --ignoz2 --numberc2 --trate2 --row2 --signh2 --group2 --scalx2 --scaly2
                    --bar2 --barstroke2 --barfill2 --opacity2  --edgey2
                    
 --help             output help information to screen

=head1 Notice
 
 1 --rankx will failure when --fredb while --group will failure without --fredb.
 2 --signs --size_st --sym_xy only work when --symbol.
 3 when use --x_scale the scale number in --scale can not less than x-ascale number, wihle scale to be a
   file the form is: posion  text in eatch line.
 4 when --samex only one x-axis value allowed, otherwise the number of x-axis and y-axis value must equal.
 5 when --bar or--bar2 the width of bar is the distance between the first two point in x-scale, so the data
   should be equidistant.
 6 --trate set is not exactly, it will be adjusted by the process.
 7 in option -ranky -rankx, the file turn start from 0 while rank from 1.
 8 option --prex --prey ordinary use 1..5, the more smaller it set, the more deatil scale will be.

=head1 Example

 perl line_diagram.pl   infile >out.svg
 read more example at /ifs2/ANIMAL/GROUP/group004/liuwenbin/SVG/example/line/output/work.sh

=cut

use strict;
use Getopt::Long;
my ($rankx,   $ranky,   $group,  $rever,   $height,  $width,   $flank_x,
    $flank_y, $size_sg, $size_d, $size_xt, $size_xs, $size_ys, $size_yt,
    $size_st, $x_mun,   $y_mun,  $linew,   $help,    $h,$size_zo,
    $x_title, $y_title, $scalx, $scaly,  $numberc, $igz,		$miceq,
    $dot,     $noline,  $color, $samex,  $fredb,   $row,	$sym_frame,
    $signh,   $sym_xy,  $signs, $symbol, $x_scale, $golden,	$scalmid,
    $valgp,   $prex,    $prey,  $frame,  $h_title, $border, $syml,
    $vice,   $valgp2, $rankx2, $ranky2,   $ankx2, $y_title2, $onedot,
    $y_mun2, $prey2,  $igz2,   $numberc2, $row2,  $signh2, $windl,
    $fredb2, $group2, $scalx2, $scaly2,   $gridx, $gridy,
    $micrx,  $micry,  $bar,    $bar2,			$micrf, $scalen,
    $barfill, $barstroke, $opacity, $barfill2, $barstroke2, $opacity2,
    $dotsel,  $nolsel,    $barsel,  $trate,    $trate2,			$size_h,
    $edgex,		$edgey,			$edgey2,    $plot, $fontfam, $italic
);
GetOptions(
           "fredb"        => \$fredb,			"size_zo:f"		=> \$size_zo,
           "samex"        => \$samex,			"valgp"       => \$valgp,
           "group:i"      => \$group,			"height:i"    => \$height,
           "linew:i"      => \$linew,			"width:i"     => \$width,
           "flank_x:i"    => \$flank_x,		"flank_y:i"   => \$flank_y,
           "size_sg:i"    => \$size_sg,		"size_xt:i"   => \$size_xt,
           "size_xs:i"    => \$size_xs,		"size_ys:i"   => \$size_ys,
           "size_yt:i"    => \$size_yt,		"size_st:i"   => \$size_st,
           "x_scale:s"    => \$x_scale,		"size_h:s"		=> \$size_h,
           "x_mun:s"      => \$x_mun,			"y_mun:s"     => \$y_mun,
           "x_title:s"    => \$x_title,		"y_title:s"   => \$y_title,
           "h_title:s"		=> \$h_title,		"border"			=> \$border,
           "rever"        => \$rever,			"ranky:s"     => \$ranky,
           "bar"          => \$bar,				"micrx"       => \$micrx,
           "micry"        => \$micry,			"micrf"				=> \$micrf,
           "miceq:s"			=> \$miceq,			"scalen"			 => \$scalen,
           "bar2"         => \$bar2,			"rankx:s"      => \$rankx,
           "sym_xy:s"     => \$sym_xy,		"signh"        => \$signh,
           "signs:s"      => \$signs,			"scalx:i"      => \$scalx,
           "gridy"        => \$gridy,			"scaly:i"      => \$scaly,
           "ignoz"        => \$igz,				"gridx"        => \$gridx,
           "dot"          => \$dot,				"color:s"      => \$color,
           "frame"        => \$frame,			"size_d:f"     => \$size_d,
           "noline"       => \$noline,		"symbol"       => \$symbol,
           "golden"       => \$golden,		"row"          => \$row,
           "prex:i"       => \$prex,			"prey:i"       => \$prey,
           "vice"         => \$vice,			"valgp2"       => \$valgp2,
           "ranky2:s"     => \$ranky2,		"rankx2:s"     => \$rankx2,
           "y_title2:s"   => \$y_title2,	"barfill:s"    => \$barfill,
           "barstroke:s"  => \$barstroke,	"opacity:i"    => \$opacity,
           "barfill2:s"   => \$barfill2,	"barstroke2:s" => \$barstroke2,
           "opacity2:i"   => \$opacity2,	"y_mun2:s"     => \$y_mun2,
           "prey2:i"      => \$prey2,			"ignoz2"       => \$igz2,
           "numberc2"     => \$numberc2,	"row2"         => \$row2,
           "signh2"       => \$signh2,		"fredb2"       => \$fredb2,
           "group2:i"     => \$group2,		"scalx2:i"     => \$scalx2,
           "scaly2:i"     => \$scaly2,		"dotsel:s"     => \$dotsel,
           "nolsel:s"     => \$nolsel,		"barsel:s"     => \$barsel,
           "trate:i"      => \$trate,		"trate2:i"     => \$trate2,
           "numberc"      => \$numberc,		"help"          => \$help,
           "sym_frame"  => \$sym_frame,	    "scalmid"		=> \$scalmid,
           "h"  => \$h,	 "edgex" => \$edgex, "italic" => \$italic,
           "edgey" => \$edgey, "edgey2" => \$edgey2, "onedot" => \$onedot,
           "plot:s" => \$plot, "syml:i"=>\$syml, "fontfam:s"=>\$fontfam,
           "windl:f"=> \$windl
);
die `pod2text $0` if ($help);
if (@ARGV < 1 || $h){
    die "Name: line_diagram.pl
Author: Wenbin Liu, liuwenbin\@genomics.org.cn
Version: 3.1,  Date: 2011-07-26
Usage: perl line_diagram.pl <infiles> [-Option] >out.svg
 <infiles>          files store Data for darwing, the data can store in more then one file
 -Option:
 1 about statistics
  -fredb            to calculate frequence distribution data to draw figure
  -group <num>      group number for frequence calculated, used when -fredb, default=50
  -windl <flo>      wind length for frequence stated, default stated acroding to -group
  -valgp            one value one group when -fredb, then -group unused
  -numberc          y-axis data not to be frequency but numberc when -fredb
  -ranky <str>      y-axis data source file0:rank1,rank2..;filei:rankj,rankk.., default 0:2
  -rankx <str>      x-axis data source, form as --ranky, default 0:1, not used when -fredb
  -samex            all the y-axis stats data use the same x-asiz data
  -row              the data store in rows instead of in ranks
  
 2 about drawing
  -dot              show the key points at the line of all group
  -noline           not show the line of all group
  -bar              to draw the bar chart instead of line chart
  -dotsel <str>     selected line show dot, '1' means the first line, unused when -dot
  -nolsel <str>     selected line not show line, '1' means the first line, unused when -noline
  -barsel <str>     selected line show bar, '1' means the first line, unused when -bar
  -barfill <str>    the bar fill colors, default set auto by the process
  -barstroke <str>  the bar stroke colors, default the same as the fill color
  -opacity <str>    the percentage of bar stroke-opacity, default=100
  -x_title <str>    the figure x-area title,default no title
  -y_title <str>    the figure y-area title,default no title
  -h                output brief help information to screen
Note: you can used -help to get detail help information\n\n";
}
foreach (@ARGV){(-f $_) || die "Error: $_ is not a file, please chack it\n";}

#========================================#
#                  MAIN                  #
#========================================#
##======        STATE      ======##
$group  ||= 50;
$group2 ||= 50;
$rankx  ||= "0:1";
$ranky  ||= "0:2";
$prex   ||= 4;
$prey   ||= 3;
$trate  ||= 95;
$trate /= 100;
$trate2 ||= 95;
$trate2 /= 100;
$fontfam ||= 'Arial';
my (@X, @Y, @X2, @Y2, @symbols);
my ($xa_min1, $xa_unit1, $xa_num1, $xa_leng1, $ya_min,  $ya_unit,  $ya_num,  $ya_leng)
  = dstat_XY( \@X,\@Y,$fredb,$valgp, $rankx, $ranky, $scalx,$scaly, $row,  $symbol, $signh,
   $x_mun,$y_mun, $numberc,$prex, $prey, $group,  $samex, $rever, $trate, $edgex, $edgey, $windl);    #sub1
my ($xa_min2, $xa_unit2, $xa_num2, $xa_leng2,$ya_min2, $ya_unit2, $ya_num2, $ya_leng2);
($vice && !$ranky2) && die "Error: when uesd --vice, --ranky2 must be set\n";
if ($vice){
    ($xa_min2, $xa_unit2, $xa_num2, $xa_leng2,$ya_min2, $ya_unit2, $ya_num2, $ya_leng2)
      = dstat_XY(\@X2,\@Y2,$fredb2, $valgp2,$rankx2,$ranky2, $scalx2, $scaly2, $row2,$symbol,
      $signh2, $x_mun,$y_mun2, $numberc2, $prex,$prey2,$group2, $samex,$rever,$trate2,$edgex,$edgey2,$windl);    #sub1
}
my ($xa_min, $xa_unit, $xa_num, $xa_leng) = ($x_mun || !$vice) ? ($xa_min1, $xa_unit1, $xa_num1, $xa_leng1)
  : mul_axis(min_max($xa_min1, $xa_min2,($xa_min1 + $xa_leng1),($xa_min2 + $xa_leng2)),$prex,0.98,$edgex);        #sub1.1.1 #sub 1.2.1
$xa_leng ||= $xa_unit * $xa_num;
if ($symbol){
    $signs && (@symbols = (-s $signs) ? (split/\n/,`less $signs`) : (split /,/, $signs));
    (@symbols != (@Y + @Y2)) && die"Error: the number of signs and y-axis values must be equal when use --symbol\n";
}
if($plot){
    my @hX = @X;
    my @hY = @Y;
    @X2 || (push @hX,@X2);
    @Y2 || (push @hY,@Y2);
    print_XY(\@hX,\@hY,$plot);
}
##======        DRAWING      ======##
$width   ||= 400;
$height  ||= 300;
$flank_x ||= $width / 4;
$flank_y ||= $height / 3;
my $pwidth  = $width + $flank_x * 2;     # Calculate the width of the paper
my $pheight = $height + $flank_y * 2;    # Calculate the height of the paper
$size_zo	||= 1;
$size_h   ||= 0.04 * $width * $size_zo;
$size_st  ||= 0.035 * $width * $size_zo;
$size_xt  ||= 0.035 * $width * $size_zo;
$size_xs  ||= 0.035 * $width * $size_zo;
$size_yt  ||= 0.06 * $height * $size_zo;
$size_ys  ||= 0.035 * $width * $size_zo;
$size_d   ||= 0.008 * $width;
$size_sg  ||= 0.008 * $width;
$linew    ||= 2;
$opacity  ||= 100;
$opacity2 ||= 100;
$miceq		||= '0.35,0.25';
($opacity, $opacity2) = ($opacity / 100, $opacity2 / 100);
my (%dotselh, %nolselh, %barselh);
if ($dotsel){ foreach (split /,/, $dotsel) { $dotselh{$_ - 1} = 1; } }
if ($nolsel){ foreach (split /,/, $nolsel) { $nolselh{$_ - 1} = 1; } }
if ($barsel){ foreach (split /,/, $barsel) { $barselh{$_ - 1} = 1; } }
##======   Creat A New drawing paper  ======##
use SVG;
my $w1 = "http://www.w3.org/2000/svg";
my $w2 = "http://www.w3.org/1999/xlink";
my $svg = SVG->new(width => $pwidth,height => $pheight, xmlns => $w1,"xmlns:xlink" => $w2);       
##=====  Draw the Line   ======##
my @colors = $color ? ((-s $color) ? (split/\s+|,/,`less $color`) : (split /,/, $color))
  : qw(crimson blue lightseagreen orange mediumpurple palegreen lightcoral dodgerblue lawngreen red olive indigogreen yellow fuchsia salmon mediumslateblue darkviolet purple sienna  black);
draw_line($xa_min,$xa_leng,$ya_min,$ya_leng, $dot,$noline,$samex,$size_d,$igz,0,\@X,\@Y,$bar,$barfill, $barstroke, $opacity,\%dotselh,\%nolselh,\%barselh,$onedot);          #sub4
$vice
  && draw_line($xa_min,$xa_leng,$ya_min2,$ya_leng2,$dot,$noline,$samex,$size_d,$igz2,$#Y + 1,\@X2,\@Y2,$bar2,$barfill2, $barstroke2, $opacity2,\%dotselh, \%nolselh, \%barselh,$onedot);     #sub4
$border && &draw_border($flank_y,$flank_x,$width,$height); #sub4+
##===== Draw the head title  ======##
my($h_x, $h_y) = ($flank_x + $width/2, $flank_y - $size_h/3);
$h_title && $svg->text('x',$h_x,'y',$h_y,'stroke','none', 'fill', 'black','-cdata', $h_title,'font-size',$size_h,'text-anchor', 'middle', 'font-family', $fontfam);
##===== Draw the Y axis  ======##
draw_yaxis($ya_min,$ya_unit, $ya_num, $ya_leng, $size_ys, $size_yt,$y_title, $frame,$gridy,$micry,$micrf,$scalen,$miceq,0,$fontfam,$italic);    #sub2
$vice && draw_yaxis($ya_min2, $ya_unit2, $ya_num2,  $ya_leng2,$size_ys, $size_yt,  $y_title2, $micrf,	0,0,$micry,$scalen,	$miceq,	1,$fontfam,$italic);    #sub2
##===== Draw the X axis  ======##
draw_xaxis($xa_min,  $xa_unit, $xa_num, $xa_leng, $size_xs, $size_xt,$x_title, $x_scale, $rever, $frame, $gridx, $micrx, $micrf, $scalen, $miceq, $scalmid,$fontfam,$italic);         #sub3
##===== Draw the symbols ======##
$syml ||= $width / 20;
my ($sym_x, $sym_y) = ($flank_x + $size_st, $flank_y + 2);
if ($sym_xy){
	my ($x, $y) = (split /,/, $sym_xy);
	($x =~ /^p/) ? ($x =~ s/^p//) : ($x= ($x - $xa_min) / $xa_leng);
	($y =~ /^p/) ? ($y =~ s/^p//) : ($y = ($y - $ya_min)  / $ya_leng);
	($sym_x, $sym_y) = ($flank_x +  $x * $width, $flank_y + (1 - $y) * $height);
}
$symbol
  && draw_symbol($size_sg,$size_st,$sym_x, $sym_y,$syml,$dot,$noline,\@symbols,$#Y,$bar,$bar2,$barfill,$barstroke,
  $opacity,  $barfill2, $barstroke2,$opacity2,  \%dotselh, \%nolselh, \%barselh, $sym_frame,$fontfam,$italic,
  $flank_x,$flank_y,$width,$height,$onedot);    #sub5
##==== Print out the Draw  ====##
print $svg->xmlify;

#========================================#
#                   SUB                  #
#========================================#
#sub00
############
sub print_XY
############
{
    my ($x,$y,$out) = @_;
    my @X = @{$x};
    my @Y = @{$y};
    open PLOT,">$out";
    foreach(0..$#Y){
        my @outx = $X[$_] ? @{$X[$_]} :  @{$X[0]};
        print PLOT "@outx\n@{$Y[$_]}\n\n";
    }
    close PLOT;
}
#========================#
#     state sub
#========================#
#sub1
#============#
sub dstat_XY
#============#
{
	my ($x,$y,$fredb,$valgp, $rankx, $ranky, $scalx,$scaly, $row,  $symbol, $signh, $x_mun, $y_mun, $numberc,
		$prex,  $prey, $group,  $samex, $rever, $rate, $edgex, $edgey, $windl) = @_;
	my ($xa_min, $xa_unit, $xa_num, $xa_leng,$ya_min, $ya_unit, $ya_num, $ya_leng);
	my (@X, @Y);
	if ($fredb){
		read_data(\@X, $ranky, $scalx, 1, $row, $symbol, $signh);
		($xa_min, $xa_unit, $xa_num, $xa_leng,$ya_min, $ya_unit, $ya_num, $ya_leng) =
		$valgp ? valgp_fredb(\@X, \@Y, $x_mun, $y_mun, $numberc, $prex, $prey, $rate, $edgex)
			: stat_fredb(\@X,\@Y,$x_mun,$y_mun,$numberc,$group, $prex, $prey,  $rate, $edgex, $windl);    #sub1.1 #sub1.2
	}else{
		read_data(\@X, $rankx, $scalx, 0, $row, $symbol, $signh);    #sub1.3
		($samex && @X != 1) && die "Error: when use --samex only one x-axis value allowed\n";
		read_data(\@Y, $ranky, $scaly, 1, $row, $symbol, $signh);    #sub1.3
		(!$samex && @X != @Y) && die"Error: when no --samex the number x-axis and y-asiz value must equaled\n";
		($xa_min, $xa_unit, $xa_num, $xa_leng,$ya_min, $ya_unit, $ya_num, $ya_leng) = 
		stat_axis(\@X, \@Y, $x_mun, $y_mun, $prex, $prey, $rate, $edgex, $edgey);    #sub1.4
	}
	$rever && (axis_rever($xa_min, $xa_leng, @X));                       #sub1.5
	@{$x} = @X;
	@{$y} = @Y;
	($xa_min, $xa_unit, $xa_num, $xa_leng, $ya_min, $ya_unit, $ya_num,$ya_leng);
}
#sub1.1
#==============#
sub valgp_fredb
#==============#
{
	#usage:($xa_min,$xa_unit,$xa_num,$xa_leng,$ya_min,$ya_unit,$ya_num,$ya_leng)=valgp_fredb(\@inX,\@outY,x_mun,y_mun,numberc,prex,prey,rate,edgex)
	my ($inX,$outY,$x_mun,$y_mun,$numberc,$prex,$prey,$rate,$edgex) = @_;
	my (@X, @Y);
	foreach my $ary (@{$inX}){
		my %sth;
		foreach (@{$ary}) { $sth{$_} ||= 0; $sth{$_}++; }
        my (@key,@value);
        @key = sort {$a<=>$b} (keys %sth);
        foreach(@key){push @value,$sth{$_};};
		push @X, [@key];
		push @Y, [@value];
	}
	if (!$numberc){
		foreach my $i (0 .. $#Y){
			my @a  = @{$Y[$i]};
			my $su = sum(@a);     #sub1.1.0
			foreach (@a) { $_ = 100 * $_ / $su; }
			$Y[$i] = [@a];
		}
	}
	@{$_[0]} = @X;
	@{$_[1]} = @Y;
	my ($xa_min, $xa_unit, $xa_num, $xa_leng) =
	$x_mun ? (split /,/, $x_mun) : mul_axis(ml_min_max(@X), $prex, 0.98, $edgex);    #sub1.1.1
	$xa_leng ||= $xa_unit * $xa_num;
	my ($ya_min, $ya_unit, $ya_num, $ya_leng) = $y_mun ? (split /,/, $y_mun)
		: mul_axis(ml_min_max(@Y), $prey, $rate);                        #sub1.1.1
	$ya_leng ||= $ya_unit * $ya_num;
	($xa_min, $xa_unit, $xa_num, $xa_leng, $ya_min, $ya_unit, $ya_num,$ya_leng);
}
#sub1.1.0
#======#
sub sum
#======#
{
	my $sum = 0;
	foreach (@_) { $sum += $_; }
	$sum;
}
#sub1.1.1
#============#
sub mul_axis
#============#
{
	#usage: (min,unit,number)=mul_axis(min,max,precision,[rate],edge)
	my ($min_x,$max_x,$prec,$rate,$edge) = @_;
	$prec ||= 4;
	$rate ||= 0.95;
    my ($u, $n) = axis_split($max_x-$min_x, $rate, $prec); #sub1.1.1.1
    my $min = $u * int($min_x / $u);
    until($min + $u * $n >= $max_x){$n++;}
	my ($leng, $mayleng) = ($u * $n,$max_x - $min);
	$edge && ($mayleng < $leng) && (($leng = $mayleng),$n--);
	($min, $u, $n, $leng);
}
#sub1.1.1.1
#=============#
sub axis_split
#=============#
{
	  #useage: axis_spli(a[,b,c]),a is the max value in the axis scale polt
    #b is the ratio of the max value to the length of the Y axis.
    #c is for precision, often use 2,4,8,16
    die"the ratio must between 0.5 to 0.98, if not you should revise your figure.\n"  if ($_[1] > 0.98 || $_[1] < 0.5);
    my ($maxv,$rate,$preci) = @_;
    $rate ||= 0.9;
    $preci ||= 2;
    $preci = 2**$preci;
    sprintf("%1.1e", $maxv) =~ /^(.)(.*)e(.*)/;
    my $mbs = $1;    # the MSB of the max value in the plot
    my $mag = $3;    # the order of magnitude of the max value in the plot.
    $mag =~ s/^\+//;
    $mag =~ s/^0*//;
    $mag  ||= 0;
    my $k = $rate / (1 - $rate) / $preci;              # the middle value used to caclutate $min_value-
                       # -you can also change preci into 2 or 1, the y-scal will become more precision
    my $min_value;     # the min value show in y axis
    foreach(2,1,0.5,0.25,0.125,0.1,0.05){
    	$min_value = $_;
    	($mbs >= $_ * $k) && last;
    }
    $min_value = $min_value * 10**$mag;
    my $value_number = int($maxv / $min_value - 1e-16) + 1;  # the number of value show in y axis
    ($min_value, $value_number);
}
#sub1.2
#==============#
sub stat_fredb
#==============#
{
	#usage:($xa_min,$xa_unit,$xa_num,$xa_leng,$ya_min,$ya_unit,$ya_num,$ya_leng)=stat_fredb(\@inX,\@outY,x_mun,y_mun,numberc,group,prex,prey,rate,edgex)
	my ($inX,$outY,$x_mun,$y_mun,$numberc,$group,$prex,$prey,$rate,$edgex,$windl) = @_;
	my ($xa_min, $xa_unit, $xa_num, $xa_leng) = $x_mun ? (split /,/, $x_mun)
		: mul_axis(ml_min_max(@{$inX}), $prex, 0.98, $edgex);    #sub1.1.1 #sub 1.2.1
	$xa_leng ||= $xa_unit * $xa_num;
	if($windl){
        $group = int($xa_leng/$windl);
        ($xa_leng % $windl) && ($group++);
    }
	my (@Y, @X);
	my $j = @{$inX};
	my ($n, $xm, $xu, $xn);
	foreach my $i (0 .. $j - 1){
		my @a = @{${$inX}[$i]};
		my $dat_sum = @a;
		if ($samex){
			foreach (@a){
				$n = int($group * ($_ - $xa_min) / $xa_leng);
				($n < 0) && ($n = 0);
				($n >= $group) && ($n = $group - 1);
				${$Y[$i]}[$n]++;
			}
		}else{
            my $mmxa = $xa_min + $xa_leng;
            foreach(@a){
      	        /[^\d\.]/ && ($_=0);
      	        ($_<$xa_min) ? ($_ = $xa_min) : (($_>$mmxa) && ($_ = $mmxa));
             }
            ($xm, $xu, $xn) = mul_axis(min_max(@a), $prex, 0.98);    #sub1.1.1
            until ($xm >= $xa_min){$xm += $xu;$xn--;}
			until ($xm + $xu * $xn <= $xa_min + $xa_leng){$xn--;}
			my $xl = $xu * $xn;
			$xu = $xl / $group;
			foreach (0 .. $group-1){
				${$X[$i]}[$_] = $_ * $xu + $xm + $xu/2;
			}
			foreach (@a){
				$n = int($group * ($_ - $xm) / $xl);
				($n < 0) && ($n = 0);
				($n >= $group) && ($n = $group - 1);
				${$Y[$i]}[$n]++;
			}
		}
		foreach (0 .. $group - 1){
			${$Y[$i]}[$_] ||= 0;
			$numberc || (${$Y[$i]}[$_] = ${$Y[$i]}[$_] * 100 / $dat_sum);
		}
	}
	if ($samex){
		my @X0;
		$xu = $xa_leng / $group;
		foreach (0 .. $group - 1){
			$X0[$_] = $_ * $xu + $xa_min + $xu/2;
		}
		@X = ([@X0]);
	}
	@{$_[0]} = @X;
	@{$_[1]} = @Y;
	my ($ya_min, $ya_unit, $ya_num) = $y_mun ? (split /,/, $y_mun)
		: mul_axis(ml_min_max(@Y), $prey, $rate);    #sub1.1.1 #sub 1.2.1
	my $ya_leng = $ya_num * $ya_unit;
	($xa_min, $xa_unit, $xa_num, $xa_leng, $ya_min, $ya_unit, $ya_num,$ya_leng);
}
#sub1.2.1
#=============#
sub ml_min_max
#=============#
{
	my ($min, $max);
	foreach my $i (@_){
		($min, $max) = $max ? min_max($min, $max, @{$i}) : min_max(@{$i});    #sub1.2.1.1
	}
	($min, $max);
}
#sub1.2.1.1
#==========#
sub min_max
#==========#
{
	my ($min, $max) = ($_[0], $_[0]);
	foreach (@_){
		($max < $_) && ($max = $_);
		($min > $_) && ($min = $_);
	}
	($min, $max);
}
#sub1.3
#============#
sub read_data
#============#
{
	#usage: read_data(\@outarray,rank,scall,weather_yvalue[0/1],row,symbol,signh)
	my ($array, $rank, $scall, $x_y, $row, $symbol, $signh) = @_;
	foreach my $i (split /;/, $rank){
		my @a = split /:|,/, $i;
		my $j = shift @a;
		foreach my $k (@a){
			my @out;
    	if ($row){
      	my $l = ($ARGV[$j]=~/\.gz$/) ? `gzip -cd $ARGV[$j] | sed -n '${k}p'` : `sed -n '${k}p' $ARGV[$j]`;
        $l =~ s/^\s+//;
        @out = split /\s+/, $l;
      }else{
      	chomp(@out = ($ARGV[$j]=~/\.gz$/) ? `gzip -cd $ARGV[$j] | awk '{print \$$k}'` : `awk '{print \$$k}' $ARGV[$j]`);
      }
      if ($signh){
      	my $sy = shift @out;
        $x_y && $symbol && (push @symbols, $sy);
      }
      if ($scall){
      	foreach (@out) { $_ ||= 0; $_ *= $scall; }
      }
      push @{$_[0]}, \@out;
    }
  }
}
#sub1.4
#============#
sub stat_axis
#============#
{
	#usage: ($xa_min,$xa_unit,$xa_num,$xa_leng,$ya_min,$ya_unit,$ya_num,$ya_leng)=stat_axis(\@X,\@Y,$x_mun,$y_mun,prex,prey,rate)
	my ($inX,$inY,$x_mun,$y_mun,$prex,$prey,$rate,$edgex,$edgey) = @_;
	my ($xa_min, $xa_unit, $xa_num, $xa_leng) = $x_mun ? (split /,/, $x_mun) : mul_axis(ml_min_max(@{$_[0]}), $prex, 0.98, $edgex);    #sub1.1.1 #sub 1.2.1
	$xa_leng ||= $xa_unit * $xa_num;
	my ($ya_min, $ya_unit, $ya_num, $ya_leng) = $y_mun ? (split /,/, $y_mun) : mul_axis(ml_min_max(@{$_[1]}), $prey, $rate, $edgey);    #sub1.1.1 #sub 1.2.1
	$ya_leng ||= $ya_unit * $ya_num;
	($xa_min, $xa_unit, $xa_num, $xa_leng, $ya_min, $ya_unit, $ya_num, $ya_leng);
}
#sub1.5
#============#
sub axis_rever
#============#
{
	my $min  = shift;
	my $leng = shift;
	foreach (@_){
		foreach (@{$_}){$_ = 2 * $min + $leng - $_;}
	}
}
##############################
#       drawing sub
##############################
#sub2
###############
sub draw_yaxis
###############
{
	#usage: draw_yaxis($ya_min,$ya_unit,$ya_num,$ya_leng,$y_rate,$size_ys,$size_yt,$y_title,$frame,$gridy,$micry,$micrf,$scalen,$miceq,$vice)
	my ($ya_min,  $ya_unit, $ya_num, $ya_leng, $size_ys, $size_yt,$y_title, $frame,   $gridy,  $micry,   $micrf,	$scalen, $miceq, $vice,$fontfam,$italic) = @_;
	my $y_rate = $height / $ya_leng;
	my ($mic1,$mic2) = split/,/,$miceq;
	my $x_edge = $vice ? $flank_x + $width : $flank_x;
	$svg->line('x1', $x_edge, 'y1', $flank_y, 'x2', $x_edge, 'y2',$pheight - $flank_y,'stroke', 'black', 'stroke-width', $linew);
	($frame && !$vice) && $svg->line('x1',$flank_x + $width,'y1',$flank_y,'x2',$flank_x + $width,'y2',$height + $flank_y,'stroke', 'black','stroke-width', $linew);
	my $add_s = $vice ? -$size_ys : $size_ys;
	my ($x_point, $y_point, $scaley) =($x_edge - $add_s * $mic1, $flank_y + $height, $ya_min);
	my $type = $vice ? 'start' : 'end';
	my $text_long = 0;
	my $m_unit = $ya_unit / 5;
	my $scal_ten = $scalen ? -1 : 1;
	for (0 .. $ya_num){
		my $scaley_long = text_long($scaley);#sub3.1.1
		($text_long < $scaley_long) && ($text_long = $scaley_long);
		$svg->text('x',$x_point,'y',$y_point + $size_ys * $mic1, 'stroke','none','fill','black','-cdata',$scaley,'font-size', $size_ys,'text-anchor', $type,'font-family', $fontfam);
		$svg->line('x1', $x_edge, 'y1', $y_point, 'x2', $x_edge - $scal_ten * $add_s * $mic1,'y2', $y_point, 'stroke', 'black', 'stroke-width', $linew);
		($frame && !$vice && $micrf)
			&& $svg->line('x1', $pwidth - $x_edge, 'y1', $y_point, 'x2', $pwidth - $x_edge + $scal_ten * $add_s * ($scalen ? $mic1 : $mic2),'y2', $y_point, 'stroke', 'black', 'stroke-width', $linew);
		if ($micry && ($_ != $ya_num)){
			my $mic_y = $y_point - $m_unit * $y_rate;
			foreach my $m (1 .. 4){
				$svg->line('x1',$x_edge,'y1',$mic_y,'x2',$x_edge - $scal_ten * $add_s * $mic2,'y2',$mic_y,'stroke','black','stroke-width', 0.75 * $linew);
				($frame && !$vice && $micrf)
					&& $svg->line('x1',$pwidth - $x_edge,'y1',$mic_y,'x2',$pwidth - $x_edge + $scal_ten * $add_s * $mic2,'y2',$mic_y,'stroke','black', 'stroke-width', 0.75 * $linew);
				$mic_y -= $m_unit * $y_rate;
			}
		}
		$gridy && $svg->line('x1',$flank_x,'y1',$y_point, 'x2',$flank_x + $width,'y2',$y_point,'stroke', 'black','stroke-width', $linew / 2);
		$y_point -= $ya_unit * $y_rate;
		$scaley += $ya_unit;
	}
	#my $g = $svg->group("transform"=>"rotate(-90,$x_point,$y_point)");
	my $rota = $vice ? 90 : -90;
	my $g = $svg->group("transform" => "rotate($rota)",'stroke'=>'none',  'fill'=>'black');
	my $gd = $golden ? 0.382 : 0.5;
	$y_title || return(1);
	my $title_length = text_long($y_title);#sub3.1.1
	my $may_size = $height / $title_length;
	($may_size < $size_yt) && ($size_yt = $may_size);
	my @alltitle;
	my @group = ('font-size',$size_yt,'text-anchor','end','font-family', $fontfam);
	my %beital;
	if($italic && ($y_title=~/\[/) && ($y_title=~/\]/)){
		foreach($y_title=~/\[(.+?)\]/g){$beital{$_}=1;}
		$y_title =~ s/^\[//;
		@alltitle = split/\[|\]/,$y_title;
	}else{
		push @alltitle,$y_title;
	}
	my $halfl = ($vice ? -1 : 1) * $size_yt * $title_length / 2;
	($x_point, $y_point) = ($x_edge - $add_s * (1 + $mic1 + $text_long),$flank_y + $gd * $height - $halfl);
	foreach my $y_title(@alltitle){
		my @group0 = @group;
		$beital{$y_title} && (push @group0,('font-style','italic'));
		my ($x_point0, $y_point0) = $vice ? ($y_point, -$x_point) : (-$y_point, $x_point); #the certer is(0,0) insure that when convert to png the title will not move
		$g->text('x', $x_point0, 'y', $y_point0,'-cdata',$y_title, @group0);
		$y_point -= ($vice ? -1 : 1) * text_long($y_title,$size_yt);
	}
}
#sub3
###############
sub draw_xaxis
###############
{
	#usage: draw_xaxis($xa_min,$xa_unit,$xa_num,$xa_leng,$size_xs,$size_xt,$x_title,$x_scale,$rever,$frame,$gridx,$micrx,$micrf,$scalen,$miceqm, $scakmid)
	my ($xa_min,  $xa_unit, $xa_num, $xa_leng, $size_xs, $size_xt,$x_title, $x_scale, $rever,  $frame, $gridx, $micrx,	$micrf, $scalen, $miceq, $scalmid,$fontfam,$italic) = @_;
	my $x_rate = $width / $xa_leng;
	my ($mc1,$mc2) = split/,/,$miceq;
	$svg->line('x1',$flank_x,'y1',$flank_y + $height,'x2',$flank_x + $width,'y2',$height + $flank_y,'stroke','black','stroke-width', $linew);
	$frame && $svg->line('x1', $flank_x, 'y1', $flank_y, 'x2', $flank_x + $width,'y2', $flank_y, 'stroke', 'black', 'stroke-width', $linew);
	my ($y_point, $scalex) = ($flank_y + $height + 1.25 * $size_xs, $xa_min - $xa_unit);
	my $scal_ten = $scalen ? -1 : 1;
	my (@scales,@scalpos);
	foreach (0 .. $xa_num){
		$scales[$_] = $xa_min + $xa_unit*$_;
		$scalpos[$_] = $flank_x + $xa_unit*$_*$x_rate;
	}
	if($x_scale){
		if(-s $x_scale){
			chomp(@scales = `awk '{print \$2}' $x_scale`);
			chomp(@scalpos = `awk '{print \$1}' $x_scale`);
			foreach(@scalpos){$_ = $flank_x + ($_ - $xa_min) * $x_rate;}
			$scalmid && (@scales == @scalpos) && (push @scalpos,($flank_x+$width));
			$xa_num = $#scalpos;
		}else{
			@scales = split /,/, $x_scale;
			(@scales < $xa_num) && die"the scale number in --scale can not less than xaxis sclae number\n";
		}
	}
	if($rever){
		$rever = -1;
		foreach(@scalpos){$_ = 2*$flank_x + $width - $_;}
	}else{
		$rever = 1;
	}
	my $scal_len = $scalpos[1]-$scalpos[0];
	foreach (0 .. $xa_num){
		$svg->line('x1',$scalpos[$_],'y1',$flank_y + $height,'x2',$scalpos[$_],'y2',$flank_y + $height + $scal_ten * $size_xs * $mc1,'stroke','black','stroke-width', $linew);
		($frame && $micrf)
			&& $svg->line('x1',$scalpos[$_],'y1',$flank_y,'x2',$scalpos[$_],'y2',$flank_y - $scal_ten * $size_xs * ($scalen ? $mc1 : $mc2),'stroke','black','stroke-width', $linew);
		$gridx && $svg->line('x1',$scalpos[$_],'y1',$flank_y + $height,'x2', $scalpos[$_], 'y2',$flank_y,'stroke','black','stroke-width', $linew / 2);
		my $x_point = $scalpos[$_];
		if($scalmid){
			($_ == $xa_num) && next;
			(-s $x_scale) && ($scal_len = $scalpos[$_+1]-$scalpos[$_]);
			$x_point += $scal_len/2;
		}
		$svg->text('x',$x_point, 'y',$y_point,'stroke','none','fill','black','-cdata',$scales[$_],  'font-size',   $size_xs,'text-anchor', 'middle', 'font-family', $fontfam);
		($_ == $xa_num) && next;
		if ($micrx){
			my $m_unit = $rever * $scal_len / 5;
			my $mic_x = $scalpos[$_] + $m_unit;
			foreach my $m (0 .. 3){
				$svg->line('x1',$mic_x,'y1',$flank_y + $height,'x2',$mic_x,'y2',$flank_y + $height + $scal_ten * $size_xs * $mc2,'stroke','black','stroke-width', 0.75 * $linew);
				($frame && $micrf)
					&& $svg->line('x1',$mic_x,'y1',$flank_y,'x2',$mic_x,'y2',$flank_y - $scal_ten * $size_xs * $mc2,'stroke','black','stroke-width', 0.75 * $linew);
				$mic_x += $m_unit;
			}
		}
	}
	$x_title || return(1);
	my $title_length = text_long($x_title);#sub3.1.1
	my $may_size = $width / $title_length;
	($may_size < $size_xt) && ($size_xt = $may_size);
	my ($xpoint, $ypoint) = ($flank_x + $width / 2 - $title_length * $size_xt/2, $flank_y + $height + (0.75 + $mc1) * $size_xs + 1.25 * $size_xt);
	write_text($xpoint,$ypoint,$x_title,$size_xt,$fontfam,$italic);#sub3.1
}
#sub 3.1
##############
sub write_text
##############
{
	my ($xpoint,$ypoint,$x_title,$size_xt,$fontfam,$italic) = @_;
	my @alltitle;
	my @group = ('font-size', $size_xt,'stroke','none','fill','black','text-anchor','start', 'font-family', $fontfam);
	my %beital;
	if($italic && ($x_title=~/\[/) && ($x_title=~/\]/)){
		foreach($x_title=~/\[(.+?)\]/g){$beital{$_}=1;}
		$x_title =~ s/^\[//;
		@alltitle = split/\[|\]/,$x_title;
	}else{
		push @alltitle,$x_title;
	}
	foreach my $x_title(@alltitle){
		my @group0 = @group;
		$beital{$x_title} && (push @group0,('font-style','italic'));
		$svg->text('x',$xpoint, 'y',$ypoint, '-cdata', $x_title,@group0);
		$xpoint += text_long($x_title,$size_xt);#sub3.1.1
	}
}
#sub3.1.1
##############
sub text_long
##############
{
	my ($text,$size,$italic) = @_;
	$italic && ($text =~ s/\[|\]//g);
	my $leng = 0;
	$leng += ($text=~s/[Labdeghnopqu_02-9]//g) * 0.555;
	$leng += ($text=~s/1//g) * 0.48;
	#$leng += ($text=~s/[Jcksvxyz]//g) * 0.5;
	$leng += ($text=~s/[=+<>]//g) * 0.58;
	$leng += ($text=~s/[r-]//g) * 0.33;
	$leng += ($text=~s/[wCDHNRU]//g) * 0.725;
	$leng += ($text=~s/m//ig) * 0.835;
	$leng += ($text=~s/W//g) * 0.945;
	$leng += ($text=~s/[ABEKPSVXY]//g) * 0.67;
	$leng += ($text=~s/[GOQ]//g) * 0.78;
	$leng += ($text=~s/[t\.,:;!\\\/\[\]]//g) * 0.275;
	$leng += ($text=~s/f//g) * 0.265;
	$leng += ($text=~s/\*//g) * 0.38;
	$leng += ($text=~s/%//g) * 0.89;
	$leng += ($text=~s/[\(\){}]//g) * 0.24;
	$leng += ($text=~s/\"//g) * 0.355;
	$leng += ($text=~s/\'//g) * 0.19;
	$text && ($leng += length($text)/2);
	$size ? $size * $leng : $leng;
}
	



#sub 4+
###############
sub draw_border
###############
{
	my ($flank_y,$flank_x,$width,$height) = @_;
	$svg->rect('x', 0, 'y', 0, 'width', 2*$flank_x + $width, 'height', $flank_y,'stroke', 'none','fill', 'white');
	$svg->rect('x', 0, 'y', $flank_y + $height, 'width', 2*$flank_x+$width, 'height', $flank_y,'stroke', 'none','fill', 'white');
	$svg->rect('x', 0, 'y', 0, 'width', $flank_x, 'height', 2*$flank_y + $height,'stroke', 'none','fill', 'white');
	$svg->rect('x', $flank_x + $width, 'y', 0, 'width', $flank_x, 'height', 2*$flank_y+$height,'stroke', 'none','fill', 'white');
}
#sub4
##############
sub draw_line
##############
{
	#usage: draw_line($xa_min,$xa_leng,$ya_min,$ya_leng,$dot,$noline,$samex,$size_d,\@X,\@Y,$n)
	 my ($xa_min,  $xa_leng, $ya_min, $ya_leng, $dot,$noline,  $samex, $size_d, $igz,$n,
		$xa, $ya, $bar,$barfill, $barstroke,$opacity, $dsh,$nsh,$bsh,$onedot) = @_;
	my ($x_rate, $y_rate) = ($width / $xa_leng, $height / $ya_leng);
	my @X = @{$xa};
	my @Y = @{$ya};
	foreach my $i (0 .. $#{$ya}){
		$i = $#Y - $i;    #reverse the lines drawing turn
		my $j = $samex ? 0 : $i;
		my @x = @{$X[$j]};
		my @y = @{$Y[$i]};
		my ($w, $h) = (($x[1] - $x[0]) * $x_rate, ($y[0] - $ya_min) * $y_rate);
		my ($x1, $y1) = ($flank_x + ($x[0] - $xa_min) * $x_rate, $flank_y + $height - $h);
		my $scl = ($barstroke || $colors[$i + $n]);
		my $fcl = ($barfill   || $colors[$i + $n]);
		if ($bar || ${$bsh}{$i + $n}){
			$svg->rect('x',$x1 - $w / 2,'y',$y1,'width',$w,'height',$h,'fill-opacity', $opacity,'stroke',$scl,'fill',$fcl);
		}else{
			($dot || ${$dsh}{$i + $n}) && make_sign($i + $n, $x1, $y1, $size_d, $onedot);    #sub4.1
		}
		foreach (1 .. $#x){
            ($y[$_] && $y[$_]=~/\d\.?/) || ($y[$_] = 0);
			$igz && (!$y[$_]) && next;
			$h = ($y[$_] - $ya_min) * $y_rate;
			my ($x2, $y2) = ($flank_x + ($x[$_] - $xa_min) * $x_rate,$flank_y + $height - $h);
			if ($bar || ${$bsh}{$i + $n}){
				$svg->rect('x',$x2 - $w / 2,'y',$y2,'width',$w,'height', $h,'stroke',$scl,'fill',$fcl,'fill-opacity', $opacity);
			}else{
				($noline || ${$nsh}{$i + $n})
					|| $svg->line('x1', $x1,'y1',$y1,'x2',$x2,'y2',$y2,'stroke',$colors[$i + $n],'fill',$colors[$i + $n],'stroke-width', $linew);
				($dot || ${$dsh}{$i + $n}) && make_sign($i + $n, $x2, $y2, $size_d, $onedot);    #sub4.1
			}
			($x1, $y1) = ($x2, $y2);
		}
	}
}
#sub4.1
##############
sub make_sign
##############
{
	# this function usde to make the signs of the different value
	#usage: make_sign(trun,cx,cy,sig_size,color_turn)
	my $sig_size = ($_[3] || 0.008 * $width);    #grobal value
	#my ($a, $b) = (($_[0] % 8), ($_[4] || int($_[0] / 8)));
	my ($a, $b) = (($_[0] % 8), int($_[0] / 8));
	my $cln   = @colors;                         #grobal value
	my $color = $colors[($a + $b) % $cln];
	$_[4] && ($a = 0);
	if ($a == 0){
		$svg->circle('cx',$_[1],'cy',$_[2],'r',$sig_size,"fill",$color,'stroke', $color);
	}elsif($a == 1){
		$svg->polygon('points',[$_[1] - $sig_size, $_[2],$_[1],$_[2] - $sig_size,$_[1] + $sig_size, $_[2],$_[1],$_[2] + $sig_size],
		'fill', $color, 'stroke', $color);
	}elsif ($a == 2){
		 $svg->polygon('points',[$_[1] - $sig_size,$_[2] + $sig_size,$_[1] + $sig_size,$_[2] + $sig_size, $_[1],$_[2] - 1.155 * $sig_size],
		'fill', $color, 'stroke', $color);
	}elsif ($a == 3){
		$svg->polygon('points',[$_[1] - $sig_size,$_[2] - $sig_size,$_[1] - $sig_size,$_[2] + $sig_size,$_[1] + $sig_size,
			$_[2] + $sig_size,$_[1] + $sig_size,$_[2] - $sig_size],'fill', $color, 'stroke', $color);
	}elsif ($a == 4){
		$svg->polygon('points',[$_[1] - $sig_size, $_[2],$_[1],$_[2] - $sig_size,$_[1] + $sig_size, $_[2],$_[1],$_[2] + $sig_size],
			'fill', $color, 'stroke', $color);
	}elsif ($a == 5){
		$svg->polygon('points',[$_[1] - $sig_size,$_[2] - $sig_size,$_[1] - $sig_size,$_[2] + $sig_size,
			$_[1] + $sig_size,$_[2] + $sig_size,$_[1] + $sig_size,$_[2] - $sig_size],'fill', 'none', 'stroke', $color);
	}elsif($a == 6){
		$svg->polygon('points',[$_[1] - $sig_size,$_[2] + $sig_size,$_[1] + $sig_size,$_[2] + $sig_size,$_[1],$_[2] - 1.155 * $sig_size],
			'fill', 'none', 'stroke', $color);
	}elsif ($a == 7){
		$svg->circle('cx',$_[1],'cy',$_[2],'r',$sig_size,"fill", 'none','stroke',$color);
	}
}

#sub5
################
sub draw_symbol
################
{
	#usage: draw_symbol($size_sg,$size_st,$sym_x,$sym_y,$leng,$dot,$noline,\@symbols)
	my ($size_sg,$size_st,$sym_x,$sym_y,$leng,$dot,$noline,$a,$n,$bar,$bar2,$barfill,$barstroke,$opacity,
		$barfill2,$barstroke2, $opacity2, $dsh,$nsh,$bsh,$sym_frame,$fontfam,$italic,
		$flank_x,$flank_y,$width) = @_;
	my @sym = @{$a};
	my $sig_long = 0;
	my $h = 2 * $size_st / 3;
	foreach (@sym){
		my $long = text_long($_,0,$italic);#sub3.1.1
		($long > $sig_long) && ($sig_long = $long);
	}
	my $sign_width = $leng + ($sig_long + 1.5) * $size_st;###
	my $max_width = $flank_x+$width - $sym_x;
	if($sym_x < $flank_x + $width && $sym_y > $flank_y && $sign_width >= $max_width){
		$size_st = ($max_width - $leng)/($sig_long+2);
		$sign_width = $leng + ($sig_long + 1.5) * $size_st;###
	}
	my $sign_height = (1.2 * @sym + 0.4) * $size_st;
	$sym_frame && $svg->rect('x',$sym_x,'y',$sym_y,'width',$sign_width,'height',$sign_height,'stroke','black','fill','none','stroke-width',$linew);
	$sym_x += 0.5 * $size_st;###
	$sym_y += 1.2 * $size_st;
	my $text_x = $sym_x + 0.5*$size_st + $leng;###
	foreach (0 .. $#sym){
		if ($bar || ${$bsh}{$_}){
			my $scl = ($barstroke || $colors[$_]);
			my $fcl = ($barfill   || $colors[$_]);
			$svg->rect('x',$sym_x,'y',$sym_y - $h,'width',$leng,'height',$h,'fill-opacity', $opacity,'stroke',$scl,'fill', $fcl);
		}else{
			($dot || ${$dsh}{$_}) && make_sign($_,$sym_x + $leng / 2,$sym_y - 0.4 * $size_st, $size_sg, $onedot);    #sub4.1
			($noline || ${$nsh}{$_})
				|| $svg->line('x1',$sym_x,'y1',$sym_y - 0.4 * $size_st,'x2',$sym_x + $leng,'y2',$sym_y - 0.4 * $size_st,'stroke',$colors[$_],'stroke-width', $linew);
		}
		write_text($text_x,$sym_y,$sym[$_],$size_st,$fontfam,$italic);#sub3.1
		$sym_y += 1.2 * $size_st;
		($_ == $n && $bar2) && ($bar = $bar2, $barfill = $barfill2, $barstroke = $barstroke2,$opacity = $opacity2);
	}
}

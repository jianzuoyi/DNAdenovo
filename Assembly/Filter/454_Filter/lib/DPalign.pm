package DPalign;
use strict;
use Exporter;
our @ISA    = qw(Exporter);
our @EXPORT = qw(global_linear local_linear global_affine local_affine);

# scoring scheme
my $m = 2;    # match
my $v = 1;    # mismatch
my $d = 3;    # gap opening
my $e = 1;    # gap extension

# B matrix: 0 -> diagonal, 1 -> up, 2 -> left
# For affine: M:0, I:1
sub global_linear {
    my $x = shift;
    my $y = shift;
    my @F;
    my @B;

    my $xl = length $x;
    my $yl = length $y;

    # initial
    $F[0][0] = 0;
    $B[0][0] = "undef";
    for ( my $i = 1 ; $i <= $xl ; $i++ ) {
        $F[$i][0] = -$i * $d;
        $B[$i][0] = 1;          # from up
    }
    for ( my $j = 1 ; $j <= $yl ; $j++ ) {
        $F[0][$j] = -$j * $d;
        $B[0][$j] = -1;         # from left
    }

    # fill
    for ( my $i = 1 ; $i <= $xl ; $i++ ) {
        for ( my $j = 1 ; $j <= $yl ; $j++ ) {
            my ( $t1, $t2, $t3 );
            if ( substr( $x, $i - 1, 1 ) eq substr( $y, $j - 1, 1 ) ) {
                $t1 = $F[ $i - 1 ][ $j - 1 ] + $m;
            }
            else {
                $t1 = $F[ $i - 1 ][ $j - 1 ] - $v;
            }
            $t2 = $F[ $i - 1 ][$j] - $d;
            $t3 = $F[$i][ $j - 1 ] - $d;

            $F[$i][$j] = $t1;
            $B[$i][$j] = 0;

            if ( $t2 > $F[$i][$j] ) {
                $F[$i][$j] = $t2;
                $B[$i][$j] = 1;
            }
            if ( $t3 > $F[$i][$j] ) {
                $F[$i][$j] = $t3;
                $B[$i][$j] = -1;
            }
        }
    }

    # traceback
    my $i   = $xl;
    my $j   = $yl;
    my $aln = traceback1( $x, $y, \@F, \@B, $i, $j );
    return $aln;
}

sub local_linear {
    my $x = shift;
    my $y = shift;
    my @F;
    my @B;

    my $xl = length $x;
    my $yl = length $y;

    # initial
    $F[0][0] = 0;
    $B[0][0] = "undef";
    for ( my $i = 1 ; $i <= $xl ; $i++ ) {
        $F[$i][0] = 0;
        $B[$i][0] = "undef";
    }
    for ( my $j = 1 ; $j <= $yl ; $j++ ) {
        $F[0][$j] = 0;
        $B[0][$j] = "undef";
    }

    # fill
    for ( my $i = 1 ; $i <= $xl ; $i++ ) {
        for ( my $j = 1 ; $j <= $yl ; $j++ ) {
            my ( $t1, $t2, $t3 );
            if ( substr( $x, $i - 1, 1 ) eq substr( $y, $j - 1, 1 ) ) {
                $t1 = $F[ $i - 1 ][ $j - 1 ] + $m;
            }
            else {
                $t1 = $F[ $i - 1 ][ $j - 1 ] - $v;
            }
            $t2 = $F[ $i - 1 ][$j] - $d;
            $t3 = $F[$i][ $j - 1 ] - $d;

            $F[$i][$j] = 0;
            $B[$i][$j] = "undef";
            if ( $t1 > $F[$i][$j] ) {
                $F[$i][$j] = $t1;
                $B[$i][$j] = 0;
            }
            if ( $t2 > $F[$i][$j] ) {
                $F[$i][$j] = $t2;
                $B[$i][$j] = 1;
            }
            if ( $t3 > $F[$i][$j] ) {
                $F[$i][$j] = $t3;
                $B[$i][$j] = -1;
            }
        }
    }

    # optimal local alignment
    my $maxF = 0;
    my $maxi = 0;
    my $maxj = 0;
    for ( my $i = 1 ; $i <= $xl ; $i++ ) {
        for ( my $j = 1 ; $j <= $yl ; $j++ ) {
            if ( $F[$i][$j] > $maxF ) {
                $maxF = $F[$i][$j];
                $maxi = $i;
                $maxj = $j;
            }
        }
    }

    # traceback
    my $aln = traceback1( $x, $y, \@F, \@B, $maxi, $maxj );
    return $aln;
}

sub traceback1 {
    my ( $x, $y, $F, $B, $maxi, $maxj ) = @_;

    my $score = $F->[$maxi][$maxj];
    my $xaln;
    my $yaln;
    my $i = $maxi;
    my $j = $maxj;
    while ( $B->[$i][$j] ne "undef" ) {
        if ( $B->[$i][$j] == 0 ) {
            $xaln .= substr( $x, $i - 1, 1 );
            $yaln .= substr( $y, $j - 1, 1 );
            --$i;
            --$j;
        }
        elsif ( $B->[$i][$j] == 1 ) {
            $xaln .= substr( $x, $i - 1, 1 );
            $yaln .= "-";
            $a = $i;
            $b = $j;
            --$i;
        }
        elsif ( $B->[$i][$j] == -1 ) {
            $xaln .= "-";
            $yaln .= substr( $y, $j - 1, 1 );
            --$j;
        }
    }
    $i++;    # alignment start on x ?
    $j++;    # alignment start on y ?
    $xaln = reverse $xaln;
    $yaln = reverse $yaln;

    my $aln;
    $aln->{score} = $score;
    $aln->{xs}    = $i;
    $aln->{xe}    = $maxi;
    $aln->{ys}    = $j;
    $aln->{ye}    = $maxj;
    $aln->{x}     = $xaln;
    $aln->{y}     = $yaln;
    return $aln;
}

sub global_affine {
    my ( $x, $y ) = @_;
    my @M;
    my @I;
    my @B;    # 3D

    my $xl = length $x;
    my $yl = length $y;

    # initial
    $M[0][0] = 0;
    $B[0][0][0] = [ "undef", "undef" ];    # [i][j][M/I] M:0 I:1
    $I[0][0]    = -1000000;
    $B[0][0][1] = [ "undef", "undef" ];

    for ( my $i = 1 ; $i <= $xl ; $i++ ) {
        $M[$i][0] = -1000000;
        $B[$i][0][0] = [ "undef", "undef" ];
        $I[$i][0] = -$d - ( $i - 1 ) * $e;
        $B[$i][0][1] = [ 1, 1 ];           # from upI
    }
    $B[1][0][1] = [ 1, 0 ];                # from upM

    for ( my $j = 1 ; $j <= $yl ; $j++ ) {
        $M[0][$j] = -1000000;
        $B[0][$j][0] = [ "undef", "undef" ];
        $I[0][$j] = -$d - ( $j - 1 ) * $e;
        $B[0][$j][1] = [ -1, 1 ];          # from leftI
    }
    $B[0][1][1] = [ -1, 0 ];               # leftM

    # fill
    for ( my $i = 1 ; $i <= $xl ; $i++ ) {
        for ( my $j = 1 ; $j <= $yl ; $j++ ) {

            # M matrix
            my ( $t1, $t2 );
            if ( substr( $x, $i - 1, 1 ) eq substr( $y, $j - 1, 1 ) ) {
                $t1 = $M[ $i - 1 ][ $j - 1 ] + $m;
                $t2 = $I[ $i - 1 ][ $j - 1 ] + $m;
            }
            else {
                $t1 = $M[ $i - 1 ][ $j - 1 ] - $v;
                $t2 = $I[ $i - 1 ][ $j - 1 ] - $v;
            }
            $M[$i][$j] = $t1;
            $B[$i][$j][0] = [ 0, 0 ];    # diagonalM

            if ( $t2 > $M[$i][$j] ) {
                $M[$i][$j] = $t2;
                $B[$i][$j][0] = [ 0, 1 ];    # diagonalI
            }

            # I matrix
            my ( $t3, $t4, $t5, $t6 );
            $t3 = $M[$i][ $j - 1 ] - $d;
            $t4 = $I[$i][ $j - 1 ] - $e;
            $t5 = $M[ $i - 1 ][$j] - $d;
            $t6 = $I[ $i - 1 ][$j] - $e;

            $I[$i][$j] = $t3;
            $B[$i][$j][1] = [ -1, 0 ];       # leftM
            if ( $t4 > $I[$i][$j] ) {
                $I[$i][$j] = $t4;
                $B[$i][$j][1] = [ -1, 1 ];    # leftI
            }
            if ( $t5 > $I[$i][$j] ) {
                $I[$i][$j] = $t5;
                $B[$i][$j][1] = [ 1, 0 ];     # upM
            }
            if ( $t6 > $I[$i][$j] ) {
                $I[$i][$j] = $t6;
                $B[$i][$j][1] = [ 1, 1 ];     # upI
            }
        }
    }

    # traceback
    my $i   = $xl;
    my $j   = $yl;
    my $aln = traceback2( $x, $y, \@M, \@I, \@B, $i, $j );
    return $aln;
}

sub local_affine {
    my ( $x, $y ) = @_;
    my @M;
    my @I;
    my @B;    # 3D
    my $xl = length $x;
    my $yl = length $y;

    # initial
    $M[0][0] = 0;
    $B[0][0][0] = [ "undef", "undef" ];    # [i][j][M/I] M:0 I:1
    $I[0][0]    = 0;
    $B[0][0][1] = [ "undef", "undef" ];

    for ( my $i = 1 ; $i <= $xl ; $i++ ) {
        $M[$i][0] = 0;
        $B[$i][0][0] = [ "undef", "undef" ];
        $I[$i][0]    = -1000000;
        $B[$i][0][1] = [ "undef", "undef" ];
    }
    for ( my $j = 1 ; $j <= $yl ; $j++ ) {
        $M[0][$j] = 0;
        $B[0][$j][0] = [ "undef", "undef" ];
        $I[0][$j]    = -1000000;
        $B[0][$j][1] = [ "undef", "undef" ];
    }

    # fill
    for ( my $i = 1 ; $i <= $xl ; $i++ ) {
        for ( my $j = 1 ; $j <= $yl ; $j++ ) {

            # M matrix
            my ( $t1, $t2 );
            if ( substr( $x, $i - 1, 1 ) eq substr( $y, $j - 1, 1 ) ) {
                $t1 = $M[ $i - 1 ][ $j - 1 ] + $m;
                $t2 = $I[ $i - 1 ][ $j - 1 ] + $m;
            }
            else {
                $t1 = $M[ $i - 1 ][ $j - 1 ] - $v;
                $t2 = $I[ $i - 1 ][ $j - 1 ] - $v;
            }
            $M[$i][$j] = 0;
            $B[$i][$j][0] = [ "undef", "undef" ];
            if ( $t1 > $M[$i][$j] ) {
                $M[$i][$j] = $t1;
                $B[$i][$j][0] = [ 0, 0 ];    # diagonalM
            }
            if ( $t2 > $M[$i][$j] ) {
                $M[$i][$j] = $t2;
                $B[$i][$j][0] = [ 0, 1 ];    # diagonalI
            }

            # I matrix
            my ( $t3, $t4, $t5, $t6 );
            $t3 = $M[$i][ $j - 1 ] - $d;
            $t4 = $I[$i][ $j - 1 ] - $e;
            $t5 = $M[ $i - 1 ][$j] - $d;
            $t6 = $I[ $i - 1 ][$j] - $e;

            $I[$i][$j] = $t3;
            $B[$i][$j][1] = [ -1, 0 ];       # leftM
            if ( $t4 > $I[$i][$j] ) {
                $I[$i][$j] = $t4;
                $B[$i][$j][1] = [ -1, 1 ];    # leftI
            }
            if ( $t5 > $I[$i][$j] ) {
                $I[$i][$j] = $t5;
                $B[$i][$j][1] = [ 1, 0 ];     # upM
            }
            if ( $t6 > $I[$i][$j] ) {
                $I[$i][$j] = $t6;
                $B[$i][$j][1] = [ 1, 1 ];     # upI
            }
        }
    }

    # optimal local alignment
    my $maxF = 0;
    my $maxi = 0;
    my $maxj = 0;
    for ( my $i = 1 ; $i <= $xl ; $i++ ) {
        for ( my $j = 1 ; $j <= $yl ; $j++ ) {
            if ( $M[$i][$j] > $maxF ) {
                $maxF = $M[$i][$j];
                $maxi = $i;
                $maxj = $j;
            }
        }
    }

    # traceback
    my $i   = $xl;
    my $j   = $yl;
    my $aln = traceback2( $x, $y, \@M, \@I, \@B, $maxi, $maxj );
    return $aln;
}

sub traceback2 {
    my ( $x, $y, $M, $I, $B, $maxi, $maxj ) = @_;

    my $xaln;
    my $yaln;
    my $i     = $maxi;
    my $j     = $maxj;
    my $score = $M->[$i][$j] >= $I->[$i][$j] ? $M->[$i][$j] : $I->[$i][$j];
    my $k = $M->[$i][$j] >= $I->[$i][$j] ? 0 : 1;    # state M or I
    while ( $B->[$i][$j][$k][0] ne "undef" ) {
        my $prek = $B->[$i][$j][$k][1];              # precursor state is M or I
        if ( $B->[$i][$j][$k][0] == 0 ) {
            $xaln .= substr( $x, $i - 1, 1 );
            $yaln .= substr( $y, $j - 1, 1 );
            --$i;
            --$j;
        }
        elsif ( $B->[$i][$j][$k][0] == 1 ) {
            $xaln .= substr( $x, $i - 1, 1 );
            $yaln .= "-";
            --$i;
        }
        elsif ( $B->[$i][$j][$k][0] == -1 ) {
            $yaln .= substr( $y, $j - 1, 1 );
            $xaln .= "-";
            --$j;
        }
        $k = $prek;
    }
    ++$i;
    ++$j;
    $xaln = reverse $xaln;
    $yaln = reverse $yaln;

    my $aln;
    $aln->{score} = $score;
    $aln->{xs}    = $i;
    $aln->{xe}    = $maxi;
    $aln->{ys}    = $j;
    $aln->{ye}    = $maxj;
    $aln->{x}     = $xaln;
    $aln->{y}     = $yaln;
    return $aln;

}
1;

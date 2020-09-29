#!/usr/bin/perl

# convert GDMA format punch file to fitted mtpl file of format produced by MDCM code
# of O. Unke

use strict;
if(@ARGV+0 != 1){ die "usage: pun2mdcm.pl <file.pun>\n"; }

my @mtpl;
my $rank=0;
my $maxrank=0;
my $natm=0;
my $i=0;
my $f=0;
my $qtot=0;

my $punch=$ARGV[0];
open(INP,"<$punch");

while(<>){
  chomp;
  my @a=split;

  if($f==7){  # rank 4 line 2
    $f=0;
    $mtpl[$i][17]=$a[3];
    $mtpl[$i][18]=$a[1];
    $mtpl[$i][24]=$a[0];
    $mtpl[$i][25]=$a[2];
  }

  if($f==6){ # rank 4 line 1
    $f=7;
    $mtpl[$i][19]=$a[4];
    $mtpl[$i][20]=$a[2];
    $mtpl[$i][21]=$a[0];
    $mtpl[$i][22]=$a[1];
    $mtpl[$i][23]=$a[3];
  }

  if($f==5){ # rank 3 line 2
    if($rank==3){ $f=0; }else{ $f=6; }
    $mtpl[$i][10]=$a[1];
    $mtpl[$i][16]=$a[0];
  }

  if($f==4){ # rank 3 line 1
    $f=5; 
    $mtpl[$i][11]=$a[4];
    $mtpl[$i][12]=$a[2];
    $mtpl[$i][13]=$a[0];
    $mtpl[$i][14]=$a[1];
    $mtpl[$i][15]=$a[3];
  }


  if($f==3){ # rank 2
    if($rank==2){ $f=0; }else{ $f=4; }
    $mtpl[$i][5]=$a[4];
    $mtpl[$i][6]=$a[2];
    $mtpl[$i][7]=$a[0];
    $mtpl[$i][8]=$a[1];
    $mtpl[$i][9]=$a[3];
  }

  if($f==2){ # rank 1
    if($rank==1){ $f=0; }else{ $f=3; }
    $mtpl[$i][2]=$a[2];
    $mtpl[$i][3]=$a[0];
    $mtpl[$i][4]=$a[1];
  }

  if($f==1){ # rank 0
    if($rank==0){ $f=0; }else{ $f=2; }
    $mtpl[$i][1]=$a[0];
    $qtot+=$a[0];
  }

  if($a[4] eq "Rank"){
    $rank=$a[5];
    if($rank > 4) {die "only implemented up to l=4, exiting\n";}
    if($rank > $maxrank){$maxrank=$rank;}
    $natm++;
    $i=$natm-1;
    $mtpl[$i][0]=$a[0]; # atom type
    $f=1;
  }

}

print "Read multipoles for $natm atoms with max. rank $maxrank\n";

open(OUT,">converted.dat");
print OUT "Qtot $qtot\n";
print OUT "converted from GDMA punch file $punch\n";
print OUT "$maxrank\n";
$i=0;
my $ii=0;
for(my $n=0; $n<=$maxrank; $n++){
    $ii=$i;
  for(my $a=0; $a<$natm; $a++){
    $i=$ii;
    for(my $t=0; $t<2*$n+1; $t++){
      my $j=$t-$n;
      my $at=$a+1;
      $i++;
      print OUT "Q$at($n,$j)     $mtpl[$a][$i] \n";
    } 
  }
}

close(OUT);



#!/usr/bin/perl

# select every nth point from cube file grid

use strict;

my $line=1;
my $natm=0;
my @title;
my @dim;
my $n1;
my $n2;
my $n3;

my $sk=3; # take every skth point

my @pots;
my @sparse;

while(<>){

  chomp;
  my @a=split;

  if($line<=$natm+6){
    $title[$line-1]=$_;
  }
  if($line==3){
    $natm=$a[0];
  }
  if($line==4){
    $n1=$a[0];
    $dim[0][0]=$a[1]*$sk;
    $dim[0][1]=$a[2];
    $dim[0][2]=$a[3];
  }
  if($line==5){
    $n2=$a[0];
    $dim[1][0]=$a[1];
    $dim[1][1]=$a[2]*$sk;
    $dim[1][2]=$a[3];
  }
  if($line==6){
    $n3=$a[0];
    $dim[2][0]=$a[1];
    $dim[2][1]=$a[2];
    $dim[2][2]=$a[3]*$sk;
  }
  if($line>6+$natm){
    my $tt;
    foreach $tt (@a){
      $pots[@pots+0]=$tt;
    }
  }

  $line++;

}

# write header for new sparse cube file
my $nn1=1+int(($n1-1)/$sk);
my $nn2=1+int(($n2-1)/$sk);
my $nn3=1+int(($n3-1)/$sk);
for(my $n=0;$n<$natm+6;$n++){
  if($n<3 || $n>5){ print "$title[$n]\n"; }elsif($n==3){
    printf " %4i   %9.6f    0.000000    0.000000\n",$nn1,$dim[0][0];
  }elsif($n==4){
    printf " %4i    0.000000   %9.6f    0.000000\n",$nn2,$dim[1][1];
  }elsif($n==5){
    printf " %4i    0.000000    0.000000   %9.6f\n",$nn3,$dim[2][2];
  }
}

my $c=0;
for(my $i=0; $i<$n1; $i++){
  for(my $j=0; $j<$n2; $j++){
    for(my $k=0; $k<$n3; $k++){
      if($c==6){
        print "\n";
        $c=0;
      }
      if($i % $sk==0 && $j % $sk==0 && $k % $sk==0){
        printf "%13.5E",$pots[$i*$n2*$n3+$j*$n3+$k];
        $c++;
      }
    }
    if($j % $sk==0 && $i % $sk==0){
      print "\n";
      $c=0;
    }
  }
}

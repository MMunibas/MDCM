#!/usr/bin/perl
#
# This script reads a list of atomic charges and an ESP CUBE file in Gaussian
# format, as well as a list of atom indices. The charges of the atoms specified
# by these indices are used to generate a background ESP, which is subtracted
# from the cube file. A new cube file is written with the listed atoms and the
# ESP from their charges removed. 
#
# The purpose is to allow fitting of part of a molecule, such as a side chain,
# while fixing certain charges at their initial values.

use strict;

if(@ARGV+0 != 4){ die "usage: remove_frag_esp_from_cube.pl <chgfile> <dens cubefile> <pot cubefile> <atomlist>\n"; }

my $cgfile=$ARGV[0];
my $dcufile=$ARGV[1];
my $cufile=$ARGV[2];
my @frz=split(/,/,$ARGV[3]);

my @q;

open(CF,"<$cgfile");

while(<CF>){
  chomp;
  my @a=split;
  if(@a+0==2){ push(@q,$a[1]); }
}

my $qfrag=0.0; #residual charge of fragment after removing requested atoms
for(my $n=0;$n<@q+0;$n++){
  $qfrag+=$q[$n];
}
for(my $n=0;$n<@frz+0;$n++){
  $qfrag-=$q[$frz[$n]-1];
}

close(CF);
open(CIN,"<$cufile");
open(COUT,">modified.pot.cube");
open(DOUT,">modified.dens.cube");

my $l=0;
my $ix=0;
my $iy=0;
my $iz=0;
my $natm=@q+0;
my $nx;
my $ny;
my $nz;
my $dx;
my $dy;
my $dz;
my @O;
my @qq;
my @xq;
my @yq;
my @zq;

while(<CIN>){
  chomp;
  my $line=$_;
  my @a=split;
  $l++;
  
  if($l==1){ 
    print COUT "$line\n";  # comment line
    print DOUT "$line\n";  # comment line
  }elsif($l==2){
    print COUT "$line\n";  # comment line
    print DOUT " Electron density\n";
  }elsif($l==3){ # no. atoms, global origin
    if($natm != $a[0]){die "Error: no. of atoms in cube does not match charges file!\n";}
    $O[0]=$a[1];
    $O[1]=$a[2];
    $O[2]=$a[3];
    printf COUT "%5i %11.6f %11.6f %11.6f %4i\n",$natm-@frz,$a[1],$a[2],$a[3],$a[4];
    printf DOUT "%5i %11.6f %11.6f %11.6f %4i\n",$natm-@frz,$a[1],$a[2],$a[3],$a[4];
  }elsif($l==4){
    $nx=$a[0];
    $dx=$a[1];
    print COUT "$line\n";
    print DOUT "$line\n";
  }elsif($l==5){
    $ny=$a[0];
    $dy=$a[2];
    print COUT "$line\n";
    print DOUT "$line\n";
  }elsif($l==6){
    $nz=$a[0];
    $dz=$a[3];
    print COUT "$line\n";
    print DOUT "$line\n";
  }elsif($l>6 && $l<=6+$natm){ #nuclear coords
    my $tf=0;
    for(my $n=0; $n<@frz+0;$n++){
      if($l-6 == $frz[$n]){$tf=1;}
    }
    if($tf==0){
      print COUT "$line\n";
      print DOUT "$line\n";
      push(@xq,$a[2]);
      push(@yq,$a[3]);
      push(@zq,$a[4]);
      push(@qq,$q[$l-7]);
    }
  }else{ #ESP grid
    for(my $n=0; $n<@a+0; $n++){
      my @r;
      $r[0]=$O[0]+$ix*$dx;
      $r[1]=$O[1]+$iy*$dy;
      $r[2]=$O[2]+$iz*$dz;
      my $V=0;
      for(my $j=0;$j<@qq+0;$j++){ #loop of charges to get ESP at point
        my @dr;
        $dr[0]=$r[0]-$xq[$j];
        $dr[1]=$r[1]-$yq[$j];
        $dr[2]=$r[2]-$zq[$j];
        my $t=sqrt($dr[0]**2+$dr[1]**2+$dr[2]**2);
        if($t<0.0001){$t=0.001;}
        $V+=$qq[$j]/$t;
      }
      printf COUT "%13.5E",$V;
      $iz++;
      if($iz>=$nz){
        $iz=0;
        $iy++;
        if($iy>=$ny){
          $iy=0;
          $ix++;
        }
      }
    }
    print COUT "\n";
  }

  
}

close(CIN);
close(COUT);

$l=0;
open(DIN,"<$dcufile");
while(<DIN>){
  chomp;
  my $line=$_;
  my @a=split;
  $l++;

  if($l>6+$natm){
    print DOUT "$line\n";
  }
  
}
close(DOUT);

printf "Residual fragment charge: %12.9f\n",$qfrag;


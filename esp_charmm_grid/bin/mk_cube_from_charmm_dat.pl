#!/usr/bin/perl

# Mike Devereux, quick script to read a list of MEP values and write into Gaussian
# cube file format

use strict;

if(@ARGV+0 != 2){ die "usage: mk_cube_from_charmm_dat.pl <cube-file> <mep-file>"; }

my $tohartree=627.503;
my $tobohr=1.88973;

my $FMEP=$ARGV[1];
my $CUB=$ARGV[0];
my ($nxpts, $nypts, $nzpts);
my ($dx, $dy, $dz);
my @origin;
my @mep; # store MEP from CHARMM

# mep file written by CHARMM contains 3-line header and list of MEP values (kcal/mol)
open(MEP,"<$FMEP");
my $l=0;
while(<MEP>){

  chomp;
  my @a=split; 
  $l++;

  if($l==1){ # num x, y, z pts 
    $nxpts=$a[0]; $nypts=$a[1]; $nzpts=$a[2];
  }
  if($l==2){ # x, y, z spacing
    $dx=$a[0]*$tobohr; $dy=$a[1]*$tobohr; $dz=$a[2]*$tobohr;
  }
  if($l==3){ # cube origin
    $origin[0]=$a[0]*$tobohr; $origin[1]=$a[1]*$tobohr; $origin[2]=$a[2]*$tobohr;
  }
  if($l>3 && @a+0==1){ #MEP point
    $mep[@mep+0]=$a[0]/$tohartree;
  }
}
close(MEP);

# cube file in Gaussian format, needed for the coordinates in the header
open(CUBE,">charmm.cube"); # output cube file
print CUBE "Cube file generated from CHARMM output $FMEP\n";

open(REF,"<$CUB");
$l=0;
my $head=6; # number of lines in header
while(<REF>){

  chomp;
  $l++;
  my $line=$_;
  my @a=split;

  if($l>1 && $l<=$head){
    print CUBE "$line\n";
  }
  if($l==3){
    $head+=$a[0];
    if(abs($a[1]-$origin[0]) > 0.001 || abs($a[2]-$origin[1]) > 0.001 || abs($a[3]-$origin[2]) > 0.001){
      die "cube file origin $a[1] $a[2] $a[3] does not match CHARMM $origin[0] $origin[1] $origin[2]\nExiting!\n";
    }
  }
  if($l==4){
    if($a[0]!=$nxpts){ die "num x pts in cube ($a[0]) does not match CHARMM ($nxpts)!\n";}
    if(abs($a[1]-$dx)>0.001){ die "dx in cube file ($a[1]) does not match CHARMM ($dx)\n";}
  }
  if($l==5){
    if($a[0]!=$nypts){ die "num y pts in cube ($a[0]) does not match CHARMM ($nypts)!\n";}
    if(abs($a[2]-$dy)>0.001){ die "dy in cube file ($a[2]) does not match CHARMM ($dy)\n";}
  }
  if($l==6){
    if($a[0]!=$nzpts){ die "num z pts in cube ($a[0]) does not match CHARMM ($nzpts)!\n";}
    if(abs($a[3]-$dz)>0.001){ die "dz in cube file ($a[3]) does not match CHARMM ($dz)\n";}
  }

}
close(REF);

# begin writing ESP from CHARMM at each grid point
my $idx=0;
for(my $ix=0;$ix<$nxpts;$ix++){
  for(my $iy=0;$iy<$nypts;$iy++){
    for(my $iz=0;$iz<$nzpts;$iz++){
      my $x=$origin[0]+$ix*$dx;
      my $y=$origin[1]+$iy*$dy;
      my $z=$origin[2]+$iz*$dz;
      my $pot=$mep[$idx++];
      printf CUBE "%13.5e",$pot;
      if($iz % 6 == 5){ print CUBE "\n"; }
    }
    print CUBE "\n";
  }
}


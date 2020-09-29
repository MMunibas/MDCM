#!/usr/bin/env python3

# script combines fitted MDCM fragment files into a single molecular file for further refinement

import sys

def usage():
    print ("Usage: python3 combine-frags.py frag1.xyz frag2.xyz ... fragN.xyz")

if len(sys.argv) < 3:
  usage()
  exit()

atob=1.8897259886

frags = [] #array of fragment file names
for i in range(1,len(sys.argv)):
  frags.append(sys.argv[i])

nfrag = len(frags)

print ('INFO: combining '+str(nfrag)+' fragments')

charges = []
RMSE = []
MAE = []
MaxAE = []
n=0
totcharge=0
for i in frags:
  with open(i,'r') as fin:
    n=n+1
    line = fin.readline().split() #number of charges
    ncharge = int(line[0])
    totcharge = totcharge + ncharge
    print ("frag"+str(n)+': '+str(ncharge)+' charges')
    line = fin.readline().split() #comment
    for j in range(ncharge):
      line = fin.readline().split() #charge coords
      if len(line) != 5:
        print ("error reading file "+i)
        exit()
      charges.append([line[0],float(line[1]),float(line[2]),float(line[3]),float(line[4])])
    line = fin.readline().split() #blank line
    line = fin.readline().split() #RMSE
    RMSE.append(float(line[1]))
    line = fin.readline().split() #MAE
    MAE.append(float(line[1]))
    line = fin.readline().split() #Max. AE
    MaxAE.append(float(line[2]))
  fin.close()

# calculate approximate fitting statistics for total molecule from fragment results
totRMSE=0.0
totMAE=0.0
totMaxAE=0.0
for i in range(len(RMSE)):
  totRMSE = totRMSE + RMSE[i]
  totMAE = totMAE + MAE[i]
  totMaxAE = max(totMaxAE,MaxAE[i])
totRMSE = totRMSE / len(RMSE)
totMAE = totMAE / len(MAE)

# write combined molecule results file
with open('combined.xyz','w') as fout:
  fout.write('%i\n' % totcharge)
  fout.write('s                      x[A]                      y[A]                      z[A]                      q[e]\n'
  )
  for i in range(len(charges)):
    fout.write('%s %25.16f %25.16f %25.16f %25.16f \n' % (charges[i][0],charges[i][1],charges[i][2],charges[i][3],
        charges[i][4]))

  fout.write('\n        RMSE %23.9e kcal/mol\n' % totRMSE)
  fout.write('         MAE %23.9e kcal/mol\n' % totMAE)
  fout.write('     max. AE %23.9e kcal/mol\n\n' % totMaxAE)
  fout.write('Coordinates in bohr\n')
  fout.write('s                   x[bohr]                   y[bohr]                   z[bohr]                      q[e]\n')
  for i in range(len(charges)):
    sign = '+'
    if charges[i][4] < 0:
      sign = '-'
    fout.write('%s %25.16f %25.16f %25.16f %25.16f \n' % (sign,charges[i][1]*atob,
        charges[i][2]*atob,charges[i][3]*atob,charges[i][4]))
    


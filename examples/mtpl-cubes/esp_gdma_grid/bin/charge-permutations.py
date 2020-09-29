#!/usr/bin/env python3

# get all permutations of positive numbers that add up to $i using $NFRAG terms

import sys

def usage():
    print ("Usage: python3 charge-permutations.py <nfrags> <ncharges>")

if len(sys.argv) < 3:
  usage()
  exit()

nfrag = int(sys.argv[1])
ncharges = int(sys.argv[2])

# array to store the combinations 
# It can contain max ncharges elements 
arr = [0] * ncharges; 

# arr - array to store the combination 
# index - next location in array 
# num - given number 
# reducedNum - reduced number  
def findCombinationsUtil(arr, index, num, 
                              reducedNum, nfrag): 
  
    # Base condition 
    if (reducedNum < 0): 
        return; 
  
    # If combination is  
    # found, print it 
    if (reducedNum == 0): 
 
      if(index == nfrag): 
        for i in range(index): 
            print(arr[i], end = " "); 
        print(""); 
      return; 
 
    if(index > nfrag): #discard solutions that require too many fragments
      return;
 
    # Find the previous number stored in arr[].  
    # It helps in maintaining increasing order 
    prev = 1 if(index == 0) else arr[index - 1]; 
  
    # note loop starts from previous  
    # number i.e. at array location 
    # index - 1 
#    for k in range(prev, num + 1): 
    for k in range(1, num + 1): 
          
        # next element of array is k 
        arr[index] = k; 
  
        # call recursively with 
        # reduced number 
        findCombinationsUtil(arr, index + 1, num,  
                                 reducedNum - k, nfrag); 


# find all combinations 
findCombinationsUtil(arr, 0, ncharges, ncharges, nfrag);


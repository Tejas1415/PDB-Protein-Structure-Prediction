# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 19:49:39 2019
@author: Tejas Krishna Reddy

Main code to produce the input and output data sequences.
"""

"""
STEP 1: Retrive the seq using biopython module         SLLPKDISQ....
STEP 2: Make it a list with real names of amino-acids  ['SER', 'LEX', 'LEX', .....] 
STEP 3: Extract CA xyz coordinates for each aminoacid and allocate to the respective aminoacid. Neglect missing amino acids.
STEP 4: Fill the missing aminoacid coordinates with average xyz of the neighbouring CA atoms.
STEP 5: Build the distance matrix and store it seperatly for each input.
STEP 6: Build an input Seq matrix.
STEP 7: Build a Deep Learning model for this problem.

"""   

""" 
THis Code is to do parsing/ extract CA coordinates, extract sequence,    \
one-hot encode the sequence, and to manually calculate the distance matrix for one particular file. 

"""


# Step 1: Sequence retreival
from Bio import SeqIO
import numpy as np
import pandas as pd

output_size = 400 # Since most of protein sequences are less than 400 amino acids, we trim to 400 - for constant length of ip/op

currentID = "4kf8.pdb"

handle = open(currentID, "rU")       #open(4kf8.pdb)
 
for record in SeqIO.parse(handle, "pdb-seqres"):
    chainA = record.seq
    print(">" + record.id + "\n" + record.seq)
    break                                               # Seq A is enough for our experiment
     
handle.close()

# Step 2: Dictionary of real names of amino acids and their corresponding alphabet representation
abbrevation = {'ALA': 'A', 'ARG': 'R', 'ASN' : 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
               'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
               'MET':'M', 'PHE':'F', 'PRO':'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
               'TYR': 'Y', 'VAL': 'V', 'ASX': 'B', 'GLX': 'Z'};
               
revAbb = dict((v,k) for k,v in abbrevation.items()) # Keys and values abv were interchanged by mistake, so corrected (reverse abbrevation)
               
chainAlist = list(chainA) #convert the string to a list of individual amino acids

chainAfullforms = list(revAbb[x] for x in  chainAlist)

# STEP 3: extracting CA atoms and assigning it to the right amino acid.
pdb = open("4kf8.pdb", "r")
database = np.array(chainAfullforms)
x_coord = list();
y_coord = list();
z_coord = list();'
j=0;
count = 0;
aminoAcids = list();
for line in pdb:
    if line[:4] == 'ATOM' and line[12:16] == " CA " and line[21:22] == "A":
        current_AA = str(line[17:20].strip())            #AA = Amino Acid
        truth = 0;                                         # Creating a Bool
        while (truth==0):
            if (current_AA == chainAfullforms[j]):       #for available atoms add their xyz coordinates
                
                aminoAcids.append(chainAfullforms[j])
                x_coord.append(line[30:38].strip())
                y_coord.append(line[38:46].strip())
                z_coord.append(line[46:54].strip())
                
                j+=1;
                truth = 1;
            else:
                aminoAcids.append(chainAfullforms[j])
                print(current_AA + "     " + chainAfullforms[j])
                x_coord.append(-999.99)                 # for missing atoms place outlier in corresponding position of the amino acid.
                y_coord.append(-999.99)
                z_coord.append(-999.99)
                
                j +=1;
                
# STEP 4: Fill the outlier values with average of previous 2 elements.
def fillNaN(list1):
    list1 = np.array(list1, dtype = float).tolist()         #Converting all elements to float in the array
    for i in range(0,len(list1)):
        if(list1[i]==-999.99):
            try:
                list1[i] = (list1[i-1] + list1[i-2]) / 2
            except:
                list1[i] = 0
                
    return list1

x_coord = fillNaN(x_coord)
y_coord = fillNaN(y_coord)
z_coord = fillNaN(z_coord)

## Fix the length of arrays to 400 x 400 since most protein sequences in PDB have length less than 400.
if len(x_coord) < output_size:        # If length of sequence is less than 400 - zero padd it to 400
    out = [0]*output_size;
    x_coord[len(x_coord):output_size] = out[len(x_coord):output_size] 
    y_coord[len(y_coord):output_size] = out[len(y_coord):output_size] 
    z_coord[len(z_coord):output_size] = out[len(z_coord):output_size] 
    aminoAcids[len(aminoAcids):output_size] = out[len(aminoAcids):output_size] 
else:                         # If length of sequence is more than 400, trim the sequence to 400
    x_coord = x_coord[:output_size]
    y_coord = y_coord[:output_size]
    z_coord = z_coord[:output_size]
    aminoAcids = aminoAcids[:output_size]
    

# Build a dataframe of all these lists.['AminoAcids', 'x_coord', 'y_coord', 'z_coord']
df = pd.DataFrame()
df['AminoAcids'] = aminoAcids
df['x_coord'] = x_coord
df['y_coord'] = y_coord
df['z_coord'] = z_coord


# Step 6: create a input distance matrix with zero padding ~400 seq length.
#  Remaining tasks - saving i/p o/p

from scipy.spatial import distance_matrix

data = np.transpose([x_coord, y_coord, z_coord])
df = pd.DataFrame(data, columns = ['x_coord', 'y_coord','z_coord'], index = aminoAcids)

df1 = pd.DataFrame(distance_matrix(df.values, df.values), index=df.index, columns = df.index)


## Now creating an input matrix.

uniqueAmino = [];             # UniqueAmino has all 22 diff amino acids possible.
for key in revAbb.keys():
    uniqueAmino.append(revAbb[key])    

ipOneHot = np.zeros([len(uniqueAmino), output_size])           

for i in range(len(uniqueAmino)):
    for j in range(output_size):
        if uniqueAmino[i] == aminoAcids[j]:
            ipOneHot[i,j] = 1
            
            ### ipOneHot is the input one hot encoded sequence. 
            ### df1 is the data frame that should be the output
            

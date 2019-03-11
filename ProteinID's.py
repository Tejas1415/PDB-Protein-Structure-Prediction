# -*- coding: utf-8 -*-
"""
A Python Script to extract all the protein Id's in RCSB PDB database.

The below link will lead to the page with all Protein ID's.
https://www.rcsb.org/pdb/rest/customReport.csv?pdbids=*&customReportColumns=structureId,structureTitle,experimentalTechnique&service=wsfile&format=csv

Created on Sat Jan 26 21:31:26 2019
@author: Tejas
"""


#Copy all info from the link to a notepad file PDBid.txt
# Read the file
proteinID = open("PDBid.txt", "r")

# Creating a list with all PDBid's in the list 
proteinList = [];
for line in proteinID:
    if line.startswith('"'):
        #print (line[1:5])
        proteinList.append(line[1:5]) #starts from 0, so 1-4

# writing the list to a file for further use.
with open('extractedPDBid.txt', 'w') as f:
    for id in proteinList:
        f.write(id)
        f.write(', ')
        
        
        
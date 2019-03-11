# -*- coding: utf-8 -*-
"""
Created on Wed Jan 30 00:15:19 2019

@author: Tejas
"""
pdb = open("4kf8.pdb", "r")

count = 0;
for line in pdb:
    if line[:4] == 'ATOM' and line[12:16] == " CA " and line[21:22] == "A":
        count += 1
print(count)


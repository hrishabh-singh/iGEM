# -*- coding: utf-8 -*-
"""
Created on Sun May 15 09:50:57 2022

@author: Hrishabh
"""

import numpy as np
from motif_finder import search # need motif_finder in the same directory

# accepts 5' to 3' sequence
# replaces any instances of "T" by "U"

def clean(seq):
    seq = seq.upper()
    mRNAseq = ""
    for base in seq:
        if base =="T":
            mRNAseq += "U"
        else:
            mRNAseq += base
    seq = mRNAseq
    return seq

def stem_loops(seq,max_nt_arc=10):
    seq = clean(seq)
    comp = {"A":["U"],"U":["A","G"],"G":["C","U"],"C":["G"]} # wobble base pairs included
    n_nt_arc = 1 # initialization of the number of nt in the arc
    possible_loops = []
    while n_nt_arc <= max_nt_arc and (len(seq)-n_nt_arc > 1): # more than 10 must be very rare
        i = 1
        while i < len(seq)-2:
            str1 = seq[:i]
            str1 = str1[::-1]
            try:
                arc = seq[i:i+n_nt_arc]
                str2 = seq[i+n_nt_arc:i+n_nt_arc+len(str1)]
                str1s = [str1[:i] for i in range(1,len(str1))]
                for string in str1s:
                    try:
                        match = search(str2[:len(string)],string,comp)
                        if match != []:
                            possible_loops.append((i,string[::-1],str2[:len(string)],n_nt_arc))
                    except:
                        pass
            except:
                pass
            i += 1
        n_nt_arc += 1
    
    return possible_loops

# perhaps some filtering and ruling out of cases required!
if __name__ == "__main__":
    sequence = "AACATGTacaataataatGGAGcatgaaCATATG"
    print(stem_loops(sequence))

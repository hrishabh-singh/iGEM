# -*- coding: utf-8 -*-
"""
Created on Mon May 16 11:15:43 2022

@author: Hrishabh
"""

import numpy as np

def clean(seq): # ensure 5' to 3'
    seq = seq.upper()
    cleanseq = ""
    for base in seq:
        if base == "T":
            cleanseq += "U"
        else:
            cleanseq += base
    return cleanseq
    
def search(seq,motif,variables_dict):
    seq = clean(seq)
    motif = clean(motif) # cleaning motif and sequence
    
    motif_indices = []
    index = 0
    for base in motif:
        if base in list(variables_dict.keys()): # getting index and type of variable in motif
            motif_indices.append((index,base))
        index += 1 
    
    def checker(possible_motif,motif_indices): 
        """checking if a sequence of same length as the motif is a motif"""
        indices = [tup[0] for tup in motif_indices] # indices of variable positions
        variable = [tup[1] for tup in motif_indices] # variables at those positions
        magicsum = 0 # some engineering, see below to know more
        
        all_indices = list(np.arange(len(possible_motif),dtype=int))
        all_constant_indices = []
        for index in all_indices:
            if index not in indices:
                all_constant_indices.append(index)
        for index in all_constant_indices:
            if possible_motif[index] == motif[index]: # checking invariable part
                pass
            else:
                magicsum -= 1 # just disrupting the order, rather arbitrary
        
        counter = 0 # for index
        for index in indices:
            if possible_motif[index] in variables_dict[variable[counter]]:
                """"checks if a base at that variable index is in the list of possible 
                bases for that particular variable"""
                magicsum +=1
            counter += 1 # increments the index for variable list

        if magicsum == len(indices): # then all conditions are met
            return True
        return False 
    
    motifs = []
    for i in range(0,len(seq)-len(motif)+1):
        possible_motif = seq[i:i+len(motif)]
        if variables_dict == {}:
            if possible_motif == motif:
                motifs.append((i,possible_motif))
        elif checker(possible_motif,motif_indices):
            motifs.append((i,possible_motif))
    
    return motifs
            
            
if __name__ == "__main__":
    sequence = "AACATGTacaataataatGGAGcatgaaCATATG" # unorganised sequence
    motif = "AAK" # unorganised motif
    variables_dict = {"N":["A","U"],
                      "K":["G","C","U"],
                      "F":["A","G"]} # Ensure uppercase here, not "T" but "U"
    # keep variable_dict variable always there, even as {} 
    
    # sequence = "AAAACUUUU"
    # motif = "AAAA"
    # variables_dict = {"A":["U"],"U":["A","G"],"G":["C","U"],"C":["G"]}
    # possible use to find hairpins; uncomment to test
    print(search(sequence,motif,variables_dict))
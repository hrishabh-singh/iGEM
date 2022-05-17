# -*- coding: utf-8 -*-
"""
Created on Fri May 13 23:00:08 2022

@author: Hrishabh
"""

dict_smiles = {"A":"OC(C(C(N1C=NC2=C(N=CN=C21)N)O3)O)C3COP(O)([O-])=O",
               "U":"OC(C(C(N1C(NC(C=C1)=O)=O)O2)O)C2COP(O)([O-])=O",
               "G":"OC(C(C(N1C=NC2=C1N=C(N)NC2=O)O3)O)C3COP(O)([O-])=O",
               "C":"OC(C(C(N1C(N=C(C=C1)N)=O)O2)O)C2COP(O)([O-])=O"}

def smile_converter(seq):
    seq = seq.upper()
    newseq = ""
    for base in seq:
        if base == "T":
            newseq += "U"
        else:
            newseq += base
    seq = newseq
    smileseq = ""
    for base in seq:
        if smileseq == "":
            smileseq += dict_smiles[base]
        else:
            indexP = -(len(dict_smiles[base])-dict_smiles[base].index("P"))
            smileseq = dict_smiles[base][:indexP+2]+smileseq+dict_smiles[base][indexP+3:]
    return smileseq

# Given sequence
if __name__ == "__main__":
    seq = "AACATGTacaataataatGGAGcatgaaCATATG"
    smiles = smile_converter(seq)
    print(smiles)
            
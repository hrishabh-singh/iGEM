# -*- coding: utf-8 -*-
"""
Created on Tue May 17 11:40:51 2022

@author: Hrishabh
"""

# an arbitrary energy score given to possible hairpins
# still does not account for the effect nt sequence present around the arc

import numpy as np
import matplotlib.pyplot as plt
from possible_hairpins import stem_loops # need possible_hairpins.py and motif_finder.py in same directory

def score(seq):
    hairpins = stem_loops(seq)
    
    scores = []
    for hairpin in hairpins:
        score_counter = 0
        counter = 0 # for string index
        for base in hairpin[1]:
            if base == "A":
                score_counter += 2
            if base == "C":
                score_counter += 3
            if base == "G":
                if hairpin[2][counter] == "C":
                    score_counter += 3
                elif hairpin[2][counter] == "U":
                    score_counter += 2
            if base == "U":
                score_counter += 2 # two bonds anyways
            counter += 1
            
        # Very arbitrary, don't know relative stability
        if hairpin[3] == 4 or hairpin[3] == 5:
            score_counter *= 10
        elif hairpin[3] == 3 or hairpin[3] == 6:
            score_counter *= 8
        elif hairpin[3] == 2 or hairpin[3] == 7:
            score_counter *= 6
        elif hairpin[3] > 7:
            score_counter *= 4
        elif hairpin[3] == 1:
            score_counter *= 2
            
        scores.append((score_counter,hairpin[0],hairpin[1],hairpin[2],hairpin[3]))
        
    return scores


# sorting a list of tuples based on one parameter

def top_n(seq,n=10): # gives top 10 by default
    scores_tup = score(seq)
    scores_sorted = sorted(scores_tup,key=lambda x:x[0],reverse=True)
    return scores_sorted[:10]

def avg_energy_score(seq):
    return np.round(sum([score[0] for score in score(seq)])/len(score(seq)),2)

def graph(scores):
    arc_nt = [score[-1] for score in scores]
    energy_scores = [score[0] for score in scores]
    plt.figure()
    plt.scatter(arc_nt,energy_scores,c="b",marker="o")
    plt.xlabel("No. of nucleotides in the arc")
    plt.ylabel("Energy Score")
    plt.grid()
    
if __name__ == "__main__":
    sequence = "AACATGTacaataataatGGAGcatgaaCATATG"
    # sequence = "GGGCUUU" # uncomment to test wobble bonding
    scores = score(sequence)
    top_10 = top_n(sequence)
    graph(scores)
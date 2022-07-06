# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 09:51:25 2022

@author: Hrishabh
"""
import numpy as np
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
from ostir.ostir import run_ostir
import time

GFP = "ATGAGTAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCAACATACGGAAAACTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCGGTTATGGTGTTCAATGCTTTGCGAGATACCCAGATCATATGAAACAGCATGACTTTTTCAAGAGTGCCATGCCTGAAGGTTATGTACAGGAAAGAACTATATTTTTCAAAGATGACGGGAACTACAAGACACGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATAGAATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTTGGACACAAATTGGAATACAACTATAACTCACACAATGTATACATCATGGCAGACAAACAAAAGAATGGAATCAAAGTTAACTTCAAAATTAGACACAACATTGAAGATGGAAGCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCCACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGAGAGACCACATGGTCCTTCTTGAGTTTGTAACAGCTGCTGGGATTACACATGGCATGGATGAACTATACAAATAA"
bbGFP = "GTTTCTTCGAATTCGCGGCCGCTTCTAGatgagtaaaggagaagaacttttcactggagttgtcccaattcttgttgaattagatggtgatgttaatgggcacaaattttctgtcagtggagagggtgaaggtgatgcaacatacggaaaacttacccttaaatttatttgcactactggaaaactacctgttccatggccaacacttgtcactactttcggttatggtgttcaatgctttgcgagatacccagatcatatgaaacagcatgactttttcaagagtgccatgcctgaaggttatgtacaggaaagaactatatttttcaaagatgacgggaactacaagacacgtgctgaagtcaagtttgaaggtgatacccttgttaatagaatcgagttaaaaggtattgattttaaagaagatggaaacattcttggacacaaattggaatacaactataactcacacaatgtatacatcatggcagacaaacaaaagaatggaatcaaagttaacttcaaaattagacacaacattgaagatggaagcgttcaactagcagaccattatcaacaaaatactccaattggcgatggccctgtccttttaccagacaaccattacctgtccacacaatctgccctttcgaaagatcccaacgaaaagagagaccacatggtccttcttgagtttgtaacagctgctgggattacacatggcatggatgaactatacaaataataaTACTAGTAGCGGCCGCTGCAGGTTTCTTC"

scar = "TACTAG"

UTRprefix = "GTTTCTTCGAATTCGCGGCCGCTTCTAGAG"

dataframe = pd.read_excel("Sequences.xlsx")

dataframe["expression"] = ""
dataframe["dG_total"] = ""

# for i,seq in enumerate(dataframe["Sequence"]):
#     Seq = UTRprefix + seq + scar + GFP[:40] # upto 40 nt of the CDS
#     pos = len(UTRprefix+seq+scar) + 1
#     results = run_ostir(Seq, name = dataframe["Symbol"][i])
#     for result in results:
#         if result["start_position"] == pos:
#             print(i+1,result["start_codon"])
#             dataframe["expression"][i] = result["expression"]
#             dataframe["dG_total"][i] = result["dG_total"]
            
# for i,seq in enumerate(dataframe["Sequence"]):
#     Seq = seq + scar + GFP[:40] # upto 40 nt of the CDS, ignoring prefix
#     pos = len(seq+scar) + 1
#     results = run_ostir(Seq, name = dataframe["Symbol"][i])
#     for result in results:
#         if result["start_position"] == pos:
#             print(i+1,result["start_codon"])
#             dataframe["expression"][i] = result["expression"]
#             dataframe["dG_total"][i] = result["dG_total"]

# tick = time.time()
# for i,seq in enumerate(dataframe["Sequence"]):
#     Seq = seq + scar + GFP # entire CDS, ignoring prefix
#     pos = len(seq+scar) + 1
#     results = run_ostir(Seq, name = dataframe["Symbol"][i])
#     for result in results:
#         if result["start_position"] == pos:
#             print(i+1,result["start_codon"])
#             dataframe["expression"][i] = result["expression"]
#             dataframe["dG_total"][i] = result["dG_total"]
# tock = time.time()
# print((tock-tick)/60)

tick = time.time()
for i,seq in enumerate(dataframe["Sequence"]):
    Seq = UTRprefix + seq + scar + GFP # entire CDS, with prefix
    pos = len(UTRprefix+seq+scar) + 1
    results = run_ostir(Seq, name = dataframe["Symbol"][i])
    for result in results:
        if result["start_position"] == pos:
            print(i+1,result["start_codon"])
            dataframe["expression"][i] = result["expression"]
            dataframe["dG_total"][i] = result["dG_total"]
tock = time.time()
print((tock-tick)/60)
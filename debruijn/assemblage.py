import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys , os
import pprint

def read_fastq (fastq):
    fastq_file = open(fastq)
    lines = iter(fastq_file.readlines())
    for line in lines :
         yield next(lines)
         next(lines)
         next(lines)


def cut_kmer(seq, taille_kmer):
  
    seq = seq.strip('\n')
    for j in range(len(seq) - taille_kmer + 1):
        yield seq[j:j+taille_kmer]

def build_kmer_dic(fastq, taille_kmer):
    
    dico = {}
    for seq in read_fastq(fastq):
        for kmer in cut_kmer(seq, taille_kmer):
            print(kmer)
            if kmer not in dico:
                dico[kmer] = 0
            dico[kmer] += 1
    return dico






if __name__ == '__main__' :
    dico = build_kmer_dic("../data/eva71_two_reads.fq", 8)
    print(dico)
    print("\n")
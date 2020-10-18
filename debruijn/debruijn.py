#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import networkx as nx
import matplotlib.pyplot as plt
from operator import itemgetter
import random
random.seed(9001)
from random import randint
import statistics

__author__ = "Sleheddine Kastalli"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Sleheddine Kastalli"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Sleheddine Kastalli"
__email__ = "slehkastalli94@gmail.com"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=21, help="K-mer size (default 21)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    return parser.parse_args()

def read_fastq (fastq):
	"""Takes a fastq_file as argument.
	Returns a sequence generator"""

	fastq_file = open(fastq)
	lines = iter(fastq_file.readlines())
	for line in lines :
		yield next(lines)
		next(lines)
		next(lines)
        
def cut_kmer(seq, l):
	'''takes as arguments a list of sequences in str and a kmer size.
	returns a k_mer generator
	'''
	seq = seq.strip('\n')
	for j in range(len(seq) - l + 1):
		yield seq[j:j+l]
        
    
	
def build_kmer_dic(fastq, l):
	"""Takes a fastq file and and a kmer size.
	Returns a dictionnary with kmer as key and its occuracy as value"""
	dico_kmer = {}
	for seq in read_fastq(fastq):
		for kmer in cut_kmer(seq, l):
			print(kmer)
			if kmer not in dico_kmer:
				dico_kmer[kmer] = 0
				dico_kmer[kmer] += 1
				return dico_kmer

def build_graph(dico_kmer):
	'''take a dictionnary of kmer and return the the
	suffix / prefix tree
	'''
	graph = nx.DiGraph()
	for kmer in dico_kmer:
		node1 = kmer[:-1]
		node2 = kmer[1:]
		nb = dico_kmer[kmer]
		graph.add_edge(node1 , node2 , weight = nb)
	return graph
        
    

	
    
    
def draw_graph(graph):
	'''
	'''
	nx.draw(graph, pos=nx.spring_layout(graph))
	plt.draw()
	plt.show()


def get_starting_nodes(graph):
	'''take a graph and return a list of entry
	nodes
	'''
	list_btw = []
	for node in graph :
		pred = list(graph.predecessors(node))
		if (not pred) :
        
			list_btw.append(node)
	return list_btw

def get_sink_nodes(graph):
	'''take a graph and return a list of exit
	nodes
	'''
	list_sink = []
	for node in graph :
		succ = list(graph.successors(node))
		if (not succ) :
			list_sink.append(node)
	return list_sink

def get_contigs(graph, list_start_node, sink_nodes):
	"""Take as argument a graph, the lists of starting and ending nodes.
	    Return a list of tuple with the contig and its size.
	"""
	contigs = []
	for input_node in list_start_node:
		for output_node in sink_nodes:
			if algorithms.has_path(graph, input_node, output_node) == True :
				path = algorithms.shortest_path(graph, input_node, output_node)
				contig = path[0]
				for i in path[1:]:
					contig = contig + i[-1]
				contigs.append((contig, len(contig)))
	return contigs


def fill(text, width=80):
	"""Split text with a line return to respect fasta format"""
	return (os.linesep.join(text[i:i+width] for i in range(0, len(text), width)))


def save_contigs(contigs, fasta):
	"""Take a list of contigs and a output file name as arguments.
	    Create a fasta file in the local depositery.
	"""
	with open(fasta_file, "w") as filout:
		for j in range(len(contigs)):
			filout.write(">contig_" +  str(number + 1) + " len=" + str(contig[1])
			+ "\n" + str(fill(contig[0])) + "\n")

def std(list_values):
	"""calculate standart deviation of a list."""
	return st.stdev(list_values)


def path_average_weight(graph, path):
	"""Take a graph and a path
	 and return the average weigth.
	"""
	all_weight = 0
	for kmer in path:
		all_weight += graph.out_degree(kmer, weight = "weight")
	average_weight = all_weight / (len(path) - 1)
	return average_weight


def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
	"""Takes a graph ,a path and remove the entry and exist nodes 
	to return a clean graph.
	"""
	entry = 1
	sink = -1
	new_graph = graph
	if delete_entry_node == True:
		entry = 0
	if delete_sink_node == True:
		sink = None
	for i in path_list:
		new_graph.remove_nodes_from(i[entry:sink])
	return new_graph


def select_best_path(graph, path_list, length_list, weight_list,
	delete_entry_node = False, delete_sink_node = False):
	"""Takes a graph, a list of paths, a list of length paths, a list of average length
		of each path and the boleans deletion or not of entry and sink nodes as
		arguments.
			Returns a graph without the unwanted paths.
	"""
	max_weight = max(weight_list)
	best_path_len = 0
	best_path_index = -1
	for li, weight in enumerate(weight_list):
		if weight > max_weight:
			max_weight = weight
			best_path_len = length_list[li]
			best_path_index = li
		elif weight == max_weight:
			if best_path_len < length_list[li]:
				best_path_len = length_list[li]
				best_path_index = li
			elif best_path_len == length_list[li]:
				best_path_index = rd.choice([best_path_index, li])
	if best_path_index == -1:
	    best_path_index = rd.randint(0, len(path_list))
	graph = remove_paths(graph, path_list[:best_path_index]
	+ path_list[best_path_index + 1:], delete_entry_node, delete_sink_node)
	return graph


def solve_bubble(graph, ancestor_node, descendant_node):
	"""Take a graph, and ancestor and descendant node as arguments.
		Return a graph without the bubble between the specified nodes.
	"""
	all_paths = list(nx.algorithms.simple_paths.all_simple_paths(graph, ancestor_node,successor_node))
	graph_bubble_path = []
	graph_bubble_len_path = []
	graph_bubble_weight = []
	for path in all_paths:
		graph_bubble_path.append(path)
		graph_bubble_len_path.append(len(path))
		graph_bubble_weight.append(path_average_weight(graph, path))
	graph_new = select_best_path(graph, graph_bubble_path, graph_bubble_len_path, graph_bubble_weight)
	return graph_new


def simplify_bubbles(graph):
	"""Take a graph as argument and return a graph without bubbles."""
	bubble_list = []
	for node in graph.nodes():
		ancestor_node = [i for i in graph.predecessors(node)]
		if len(ancestor_node) > 1:

			ancestor = nx.lowest_common_ancestor(graph, ancestor_node[0], ancestor_node[1])
			bubble_nodes.append([ancestor, node])
	for node_couples in bubble_nodes:
	 	graph = solve_bubble(graph, node_couples[0], node_couples[1])
	return graph


def solve_entry_tips(graph, entry):
	"""Takes a graph and all entry nodes.
	  Returns: A graph with no entry tips

	"""
	Pass


def solve_out_tips(graph, nodes_out):
	"""Takes a directed graph and all exist nodes 
	  Returns: A graph with no outtips
	"""
	Pass
#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments

    args = get_arguments()


    PARSER = argparse.ArgumentParser()

    PARSER.add_argument("fasta_file", help="the fasta file", type=str)
    PARSER.add_argument("len_kmer", help="the length of the kmers", type=int)

    ARGS = PARSER.parse_args()

    FASTA_FILE = ARGS.fasta_file
    LEN_KMER = ARGS.len_kmer

    dico = build_kmer_dic(FASTA_FILE, LEN_KMER)
    print(dico)
    print("\n")
    G = build_graph(dico_kmer)
    draw_graph(G)
    starting_nodes = get_starting_nodes(G)
    print(starting_nodes)
    sink_nodes = get_sink_nodes(G)
    print(sink_nodes)
    G = simplify_bubbles(G)
    
    print(G)
    contigs = get_contigs(G,get_starting_nodes(G),get_sink_nodes(G))
    save_contigs(contigs,args.output_file)



    
if __name__ == '__main__':
	main()
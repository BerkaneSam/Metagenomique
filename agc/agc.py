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

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
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
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication  (default 100)")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()

def read_fasta(amplicon_file, minseqlen):
    with gzip.open(amplicon_file, 'rt')as filin:
        seq = ""
        for line in filin:
            if not line.startswith(">"):
                seq += line.strip()
            else:
                if len(seq) >= minseqlen:
                    yield seq
                seq = ""
        yield seq
                


def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    seq_gen = list(read_fasta(amplicon_file, minseqlen))
    seq_count = [(seq, seq_gen.count(seq)) for seq in seq_gen]
    seq_count = list(get_unique(seq_count))
    seq_count.sort(key=lambda item: item[1], reverse=True)
    for seq in seq_count:
        if seq[1] >= mincount:
            yield list(seq)




def get_unique(ids):
    return {}.fromkeys(ids).keys()


def common(lst1, lst2): 
    return list(set(lst1) & set(lst2))


def get_chunks(sequence, chunk_size):
    """"""
    len_seq = len(sequence)
    if len_seq < chunk_size * 4:
        raise ValueError("Sequence length ({}) is too short to be splitted in 4"
                         " chunk of size {}".format(len_seq, chunk_size))
    return [sequence[i:i+chunk_size] 
              for i in range(0, len_seq, chunk_size) 
                if i+chunk_size <= len_seq - 1]


def cut_kmer(sequence, kmer_size):
    """Cut sequence into kmers"""
    for i in range(0, len(sequence) - kmer_size + 1):
        yield sequence[i:i+kmer_size]

def get_identity(alignment_list):
    """Prend en une liste de séquences alignées au format ["SE-QUENCE1", "SE-QUENCE2"]
    Retourne le pourcentage d'identite entre les deux."""
    id_nu = 0
    for i in range(len(alignment_list[0])):
        if alignment_list[0][i] == alignment_list[1][i]:
            id_nu += 1
    return round(100.0 * id_nu / len(alignment_list[0]), 2)


def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    kmers = cut_kmer(sequence, kmer_size)
    for kmer in kmers:
        if kmer not in kmer_dict:
            kmer_dict[kmer] = [id_seq]
        else:
            kmer_dict[kmer].append(id_seq)
    return kmer_dict


def search_mates(kmer_dict, sequence, kmer_size):
    kmers = cut_kmer(sequence, kmer_size)
    id_list = []
    for kmer in kmers:
        if kmer in kmer_dict:
            id_list += kmer_dict[kmer]
    result = Counter(id_list).most_common(2)
    id_parents = [id[0] for id in result ]
    return id_parents


def detect_chimera(perc_identity_matrix):
    std_results = []
    par1 = False
    par2 = False
    for seg_id in perc_identity_matrix:
        std_results.append(statistics.stdev(seg_id))
        if seg_id[0] > seg_id[1]:
            par1 = True
        if seg_id[1] > seg_id[0]:
            par2 = True
    std_mean = statistics.mean(std_results)
    if std_mean > 5:
        if par1 == True and par2 == True:
            return True
    return False
        





def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    pass

def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    otu = []
    mother_intermediate = []
    seqs = list(dereplication_fulllength(amplicon_file, minseqlen, mincount))
    for i in range(len(seqs)):
        if seqs[i] in mother_intermediate:
            continue 
        intermediate = [seqs[i]]
        for j in range(i + 1, len(seqs)):
            seq_align = nw.global_align(seqs[i][0], seqs[j][0], gap_open=-1, gap_extend=-1, 
            matrix=os.path.abspath(os.path.join(os.path.dirname(__file__),"MATCH")))
            id = get_identity(seq_align)
            if id > 97:
                intermediate.append(seqs[j])
        otu.append(intermediate[0])
        mother_intermediate += intermediate
    return otu


def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def write_OTU(OTU_list, output_file):
    count = 1
    with open(output_file, 'w')as filout:
        for otu in OTU_list:
            filout.write(f">OTU_{count} occurrence:{otu[1]}\n")
            seq = fill(otu[0])
            filout.write(f"{seq}\n")
            count += 1


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # Votre programme ici


if __name__ == '__main__':
    main()

# -*- coding: utf-8 -*-

import pandas as pd
from Bio import SeqIO

def read_aligned_files(path, name):
    '''
    Read aligned files
    :param path: path of the file
    :param name: name of the file
    :return: return a Dataframe with the file readed
    '''
    path += name
    temp = []
    with open(path) as fasta_file:
        for seqRecord in SeqIO.parse(fasta_file, "fasta"):
            temp.append(seqRecord)
    return temp

def seq_to_df(sourceList):
    """
    Save on dataframes files of type Bio.Seq
    :param sourceList: Bio.Seq file
    :return: A dataframe
    """
    baseDf = pd.DataFrame(columns=["Name", "Seq", "Len"])

    for itr in sourceList:
        baseDf.loc[len(baseDf)] = [itr.name, str(itr.seq), len(itr.seq)]

    return baseDf

def readPDBs(dataList, option = "node"):
    """
        Read my PDBs files
        dataList is the list of pdbs names
        option default or equals to 'node', read nodes
        option = 'edges', read edges files
    """
    src = "_nodes.txt"
    if option == "edges":
        src = "_edges.txt"
        
    nsFiles = {}
    for ii in dataList:
        nsFiles[ii.upper()] = pd.read_csv("read/PDBs RIN/"+ii+src, sep="\t", low_memory=False)
    
    return nsFiles
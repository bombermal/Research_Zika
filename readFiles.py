# -*- coding: utf-8 -*-

from Bio import SeqIO
import pandas as pd

def fileRead(path, type):
    """
       Read files and return a list with all values
       Parameters
       ----------
       path: string
           string with the path of the file
        ext: string
            string with the file extension
    """
    temp = []
    with open(path) as fasta_file:
        for seqRecord in SeqIO.parse(fasta_file, type):
            temp.append(seqRecord)
    return temp

def genomeClustalToDf(columns, sourceList, organism):
    """
        Make a dataFrame
        
        Parameters
        ----------
        souceList: list
            Data to populate DF
        columns: list
            List with columns names
        
        Return
        ------
        pandas.DataFrame
    """
    
    
    baseDf = pd.DataFrame(columns=columns)
    
    for itr in sourceList:
        tam = len(baseDf)
        baseDf.loc[tam] = [organism, itr.description, str(itr.seq), len(str(itr.seq))]
    return baseDf


def genomeFastaToDf(columns, sourceList):
    """
        Make a dataFrame
        
        Parameters
        ----------
        souceList: list
            Data to populate DF
        columns: list
            List with columns names
        
        Return
        ------
        pandas.DataFrame
    """
    
    
    baseDf = pd.DataFrame(columns=columns)
    
    for itr in sourceList:
        tam = len(baseDf)
        description = itr.description.split("|")
        baseDf.loc[tam] = [description[0].replace(" ", ""), description[2], description[3]]
    return baseDf

def readPDBs(dataList, option = "node", cut = -3.7):
    """
        Read my PDBs files
        dataList is the list of pdbs names
        option default or equals to 'node', read nodes
        option = 'edges', read edges files
        option = discotop, read discotop files
    """
    src = "_nodes.txt"
    if option == "edges":
        src = "_edges.txt"
    if option == "discotop":
        src = "result_"
        
    nsFiles = {}
    for ii in dataList:
        if option == "discotop":
            nsFiles[ii.upper()] = pd.read_csv("read/PDBs/"+src+str(abs(cut))+"_"+ii.lower()+".csv", skiprows=1, sep=",", low_memory=False)
        else:
            nsFiles[ii.upper()] = pd.read_csv("read/PDBs/"+ii+src, sep="\t", low_memory=False)
    
    return nsFiles
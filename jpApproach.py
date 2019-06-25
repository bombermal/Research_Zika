# -*- coding: utf-8 -*-
##Imports
import pandas as pd
import readFiles as rd
import processFunctions as fc
import workFlow as wf
import graph as gr

#%% Prepare work.
%%time
#Read files
filesPaths = { "Zika" : "read/Zika Full 454 29-10-19.fa", "Dengue" : "read/Dengue Full 4827 29-10-19.aln"}
file = rd.fileRead(filesPaths["Zika"], "fasta")

#Convert file to a DataFrame
baseDf = rd.genomeClustalToDf(["Organism", "SeqId", "Seq", "Len"], file, "Zika")

#Create Dataframes
sequencesTable = fc.createDfForCount(baseDf)
snpsCountDf = fc.makeDfForPolimorphismsCount(sequencesTable)

#Select poliorphic positions
snpsDf = snpsCountDf[snpsCountDf["Cons/Poli"] == "P"]
poliPositionsList = snpsDf.index.tolist() 

# Read PDB's
pdbsNames = [ "5GS6", "5K6K", "5IY3", "5JMT", "5TMH"]
"""
    rd.readPDBs(fileNames, option, discotop value)
    fileName = list with pdbs names
    option = 'node', 'edges', 'discotop'
    discotop value = -3.7 or -7.7
"""
pdbsNodesFiles = rd.readPDBs(pdbsNames)
pdbsEdgesFiles = rd.readPDBs(pdbsNames, "edges")

#Split das amostras e PDB
    
samplesList = fc.splitSmapleDataFrame(baseDf) 
pdbsList = fc.splitPdbsDict(pdbsNodesFiles)
#%% Read Hit Table

def readHitFiles(pdbsNames):
    """
        Read HitTable files from Blast and save all in dicts
    """
    temp = {}
    #read files
    for key in pdbsNames:
        cols = "query acc.ver,subject acc.ver,% identity,alignment length,mismatches,gap opens,q. start,q. end,s. start,s. end,evalue,bit score,% positives".split(",")
        aux = pd.read_csv("pdbs fasta/"+key+".csv", header = None, names = cols)
        temp[key] = aux
        
    return temp

def filterHitTable(hitTableDF):
    """
        Remove the duplicated values from dictionarys
    """
    temp = {}
    for key, value in hitTableDF.items():
        print()
    

hitTableDF = readHitFiles(pdbsNames)

#filteredHitTables = filterHitTable(hitTableDF)
#%% Separar sequencias
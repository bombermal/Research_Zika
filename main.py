# -*- coding: utf-8 -*-
##Imports
import pandas as pd
import readFiles as rd
import processFunctions as fc
import workFlow as wf
import graph as gr

#%% Prepare work.
#%%time
#Read files
filesPaths = { "Zika" : "read/Zika Full 454 29-10-19.aln", "Dengue" : "read/Dengue Full 4827 29-10-19.aln"}
file = rd.fileRead(filesPaths["Zika"], "clustal")

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

#%% Alinhamento
#%%time

def readOrWork(option="Read"):
    if option == "Work":
        obj = wf.workFlow(pdbsList, samplesList)
        obj.getResults().to_csv("read/alignedDataFrame.csv", sep="\t", index=False)
    #Lê
    alignedDf = pd.read_csv("read/alignedDataFrame.csv", sep="\t")

    return alignedDf

alignedDf = readOrWork()

#%%
#%%time
minMaxPoitionsPerPdb = fc.listMinMaxTuple(pdbsNames, alignedDf)
poliPositionsDict = fc.convertPositionStringToInt(poliPositionsList, minMaxPoitionsPerPdb)    
aminoacidsOccurrencesDict = fc.totalAminoCountingPerPositions(poliPositionsDict, alignedDf)

#%% Plot cobertura
"""
Plot cover distribution
"""

condition = True
gr.plotCoverRange(alignedDf, condition)

#Show PDBs degrees distribution

for name in pdbsNames:
    """
     True = Save figures
     False = Don't save
    """
    title = "Overview of degrees distribution: "
    gr.plotDegreeDistribution(title, name, pdbsNodesFiles[name], condition)
   
gr.plotPolimorphismsDistribution(aminoacidsOccurrencesDict, pdbsNodesFiles, condition)

#%% New Approach

#Print PDBS sequences
#For blast
#fc.printSequencesForDebug(pdbsList)

def saveCSV(path, df):
    df.to_csv(path, sep="\t", index=False)      

#saveCSV("Df Polimorfismos anotados.csv", snpsDf)
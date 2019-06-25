# -*- coding: utf-8 -*-
import pandas as pd
from Bio.SeqUtils import seq3, seq1 as seqX
from collections import Counter
import node as nd

def createDfForCount(baseDf):
    """
        Convert sequences in a table, with aminoacids aligned
    """
    columns = ["SeqId"]+[str(x)+"POS" for x in range(baseDf.loc[0, "Len"])]
    aux = []
    for idx, row in baseDf.iterrows():
        aux.append([row.SeqId]+list(row.Seq))  
    
    return pd.DataFrame(aux, columns=columns)

def returnPercent(totalOfSamples, percent):
    return ((totalOfSamples*percent)/100, (totalOfSamples*abs(percent-100)/100))

def returnValuesToSNPsCountDf(val):
    """
        Count SNPs frequency by column
        
        return a list for dataFrame merge
    """
    snps = []
    frequency = []
    xFrequency = ""
    
    for amino, freq in val.iteritems():
            if amino == "X":
                xFrequency = amino+":"+str(freq)
            else:
                snps.append(amino)
                frequency.append(freq)
    aux = []
    
    for a, f in zip(snps, frequency):
        aux.append(a+":"+str(f))
        
    limit = returnPercent(454, 1) # 454 amostras/ acima de 1% de mudanças são vistas como polimorfismos
    conservation = "P" if frequency[0] < limit[1] else "C"
    
    return [",".join(snps), conservation , xFrequency, ",".join(aux)]

def makeDfForPolimorphismsCount(preContabilizationDf):
    """
        Create a Df with SNPs counted
    """
    aux = []
    columns = preContabilizationDf.columns[1:]
    for col in columns:
        temp = preContabilizationDf[col].value_counts()
        aux.append(returnValuesToSNPsCountDf(temp))
        
    return pd.DataFrame(aux, index=columns, columns=["SNPs", "Cons/Poli", "X", "Common/Frequency"])

def splitSmapleDataFrame(df):
    """
        Split the DataFrame in a list of series (rows)
    """
    temp = []
    for idx, row in df.iterrows():
        temp.append(row)
    
    return temp

def splitPdbsDict(pDict):
    """
        Divide o Df em dicionários, onde key = PdbId e value = Df
    """
    temp = {}
    for key in pDict:
        aux = {}
        for chain in pDict[key].Chain.unique():
            aux[chain] = makeNodes(pDict[key], True)
        temp[key] = aux
    return temp

def makeNodes(seq, pdb=False):
    """
        Make list of nodes from a string
        True if the file was a PDB
        False if a Sample
    """
    temp = []
    if pdb:
        for idx, amino in seq.iterrows():
            nodeId = seq.loc[idx, "NodeId"]
            temp.append(nd.node(seqX(seq.loc[idx, "Residue"]), seq.loc[idx, "Degree"], nodeId, nodeId.split(":")[1], -2))
    else:
        for idx, amino in enumerate(seq):
            nodeId = "@:"+str(idx)+":_:"+amino
            temp.append(nd.node(amino, -10, nodeId, -2 , idx))
            
    return temp 

def listMinMaxTuple(pdbsNames, alignedDf):
    """
        Returna a dict with tuples of min and max values of position encountered in alignedDf
    """
    minMaxPoitionsPerPdb = {}
    for key in pdbsNames:
        temp = alignedDf[alignedDf.PdbId == key]
        for idx, seq in temp.iterrows():
            aux = [int(x) for x in seq["SampleSeqPos"].split(",")]
            minMaxPoitionsPerPdb[key] = (min(aux), max(aux))
    
    return minMaxPoitionsPerPdb

def convertPositionStringToInt(poliPositionsList, minMaxPoitionsPerPdb):
    """
        Return a dictionary with a list of selected polimorphic positions per pdb
    """
    poliPositionsList = [int(x[:-3]) for x in poliPositionsList]
    aux = {}
    for key, value in minMaxPoitionsPerPdb.items():
        temp = []
        mini = value[0]
        maxi = value[1]
        for ii in poliPositionsList:
            if ii >= mini and ii <= maxi:
                temp.append(ii)
                
        aux[key] = temp
        
    return aux

def rowToAminoList(row, position):
    """
        Search the received row for positions wanted.
       Return a list with the aminoacids on this positions 
    """
    temp = []
    for aminoId, seqPos in zip(row.SeqAminoId.split(","), row.SampleSeqPos.split(",")):
        if int(seqPos) == position:
            temp.append(aminoId)
            
    return temp

def totalAminoCountingPerPositions(poliPositionsDict, alignedDf):
    """
        Create a dictionary with all aminoacids found in the polimorphic positions
    """
    finalTemp = {}
    for key, value in poliPositionsDict.items():
        tempDf = alignedDf[alignedDf.PdbId == key]
        temp = []
        for position in value:
            for idx, row in tempDf.iterrows():
                temp += rowToAminoList(row, position )
                
        finalTemp[key] = temp
        
    return finalTemp

def countingOcurrences(aminoacidsOccurrencesDict):
    """
        Use Counter Library to count the frequency of each aminoacid
    """
    temp = {}
    for key, value in aminoacidsOccurrencesDict.items():
        temp[key] = Counter(value)
    
    return temp
    
def addDegreeToOccurrencesCount(countDict, pdbsNodesFiles):
    """
        Create a list with degrees frequency for each PBD
        
    """
    aux = {}
    for key in countDict:
        temp = {}
        for subKey, value in countDict[key].items():
            temp[subKey] = (value, pdbsNodesFiles[key][pdbsNodesFiles[key].NodeId == subKey].Degree.values[0])
        
        aux[key] = temp
    return aux

def writeSequencesForDebug(pdbsList):
    for key in pdbsList:
        #for chain in pdbsList[key]:
        chain = "A"    
        file = open("pdbs fasta/"+key+".fa", "a")
        file.write(">{}:{}\n{}\n".format(key, chain, "".join([ii.amino for ii in pdbsList[key][chain]])))
        file.close()
            
def printSequencesForDebug(pdbsList):
    for key in pdbsList:
        for chain in pdbsList[key]:
            print(">{}:{}\n{}\n".format(key, chain, "".join([ii.amino for ii in pdbsList[key][chain]])))
            
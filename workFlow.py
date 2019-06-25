# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 01:20:57 2018

@author: Ivan Alisson
"""

import pandas as pd
from tqdm import tqdm
from multiprocessing.dummy import Pool as ThreadPool
import multiprocessing as mp
import processFunctions as fc
import smithWaterman as sw

class workFlow:
    pdbsDict = None
    samplesList = None
    results = None
    columns = ["PdbId", "SampleId", "Chain", "SeqAminoId", "SampleSeqPos", "Cover"]
    
    def __init__(self, pdbsDict, samplesList):
        self.pdbsDict = pdbsDict
        self.samplesList = samplesList
        
        self.makeAllWork()
        
    def singleStep(self, oneSample):
        tempDf = pd.DataFrame(columns=self.columns)
        seqSample = fc.makeNodes(oneSample.Seq, False)
        for key in self.pdbsDict:
            for chain in ["A"]:#self.pdbsDict[key]:
                seqPdb = self.pdbsDict[key][chain]
                obj = sw.smithWaterman()
                self.seqIds, self.seqPos, cut = obj.constructor(2, -1, -1, seqSample, seqPdb, False, False)
                #self.seqIds, self.seqPos, cut = al.water(seqSample, seqPdb)
                
                auxIds = ",".join(self.seqIds)
                auxPos = ",".join(self.seqPos)
                
                tam = len(tempDf)
                tempDf.loc[tam] = key, oneSample.SeqId, chain, auxIds, auxPos, cut
        
        return tempDf
    
    def makeAllWork(self):
        
        poolSize = mp.cpu_count()
        pool = ThreadPool(poolSize)
        results = pd.DataFrame(columns=self.columns)
        
        for ii in tqdm(pool.imap_unordered(self.singleStep, self.samplesList), total=len(self.samplesList)):
            results = pd.concat([results, ii])
        
        #for ii in tqdm(self.samplesList):
        #    results = pd.concat([results, self.singleStep(ii)])
        
        pool.close()
        pool.join()
        
        results.reset_index(inplace=True, drop=True)
        
        self.results = results
        
    def getResults(self):
        return self.results
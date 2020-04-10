# -*- coding: utf-8 -*-
"""
Created on Tue Dec  5 17:37:26 2017

@author: Ivan Alisson Cavalcante Nunes de Lima
"""

def make_node(row):
    return 1

class node:
    """
        Class tha stores each string in a node, with a character, a index before alignment and a index after the alignment
    """
    amino = None
    degree = None
    aminoId = None
    pdbPos = None
    seqPos = None
    
    def __init__(self, amino, degree, aminoId, pdbPos, seqPos):
        self.amino = amino
        self.degree = degree
        self.aminoId = aminoId
        self.pdbPos = pdbPos
        self.seqPos = seqPos

    def getAll(self):
        return self.amino, self.degree, self.aminoId, self.pdbPos, self.seqPos
    
    def setAll(self, amino, degree, aminoId, pdbPos, seqPos):
        self.amino = amino
        self.degree = degree
        self.aminoId = aminoId
        self.pdbPos = pdbPos
        self.seqPos = seqPos

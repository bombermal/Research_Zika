# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import processFunctions as fc

def plotDegreeDistribution(title, pdbName, pdb, save=False):
    fig = plt.figure()
    pdb.Degree.plot(kind="hist", title=title+pdbName, colormap="tab20")
    if save:
        fig.savefig("figure/Degrees_distribution_"+pdbName)
        
def plotCoverRange(df, save=False):
    fig = plt.figure()
    df.Cover.plot(kind="hist", title="Overview identities", colormap="tab20")
    if save:
        fig.savefig("figure/Cover distribution") 
        
def plotPolimorphismsDistribution(aminoacidsOccurrencesDict, pdbsNodesFiles, save=False):
    temp = fc.addDegreeToOccurrencesCount(fc.countingOcurrences(aminoacidsOccurrencesDict), pdbsNodesFiles)
    for key in temp:
        aux = []
        for subkey, value in temp[key].items():
            aux += [value[1] for _ in range(value[0])]
            
        fig, ax = plt.subplots(1,1)
        ax.set_title("Polimorphism degree distribution "+key)
        ax.hist(aux)
        if save:
            fig.savefig("figure/Polimorphism_Distribution_"+key)
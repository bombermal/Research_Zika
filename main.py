# Imports
import pandas as pd
from IPython import get_ipython

# My imports
import Codes.read as rd
import Codes.function as fc
import Codes.graph as gf
#%% Reading files

#Leio o arquivo alinhado
path = "read/"
name = "Zika Full 454 29-10-19.fasta"
raw_aln_zika = rd.read_aligned_files(path, name)
#Converto alinhado para um DF
raw_df_aln = rd.seq_to_df(raw_aln_zika)

#%%


def read_Or_create(raw_df_aln, read=True):
    if read:
        #lê salvo
        counted_df = pd.read_csv("Saved/1_Counting_df_zika_454_29-10-19.csv")
    else:
        #Df vazio
        counted_df = pd.DataFrame()
        #função que faz a transposta de Seq e conta as ocorrencias
        fc.transpose_seq_and_count(raw_df_aln.Seq, counted_df)
        #Limpa os Nan e corrige a Pos
        counted_df = counted_df.fillna(0).reset_index().rename(columns={"index" : "Pos"})
        counted_df.Pos = counted_df.Pos.apply(lambda x: x+1)
        #Salva o trabalho
        counted_df.to_csv("Saved/1_Counting_df_zika_454_29-10-19.csv", index=False)
    
    return counted_df
        
counted_df= read_Or_create(raw_df_aln)
#%% Ploting

conditions = [.5,1,2]
list_of_filtered_dfs = []

for ii in conditions:
  list_of_filtered_dfs.append(fc.filter_criteria(counted_df, len(raw_aln_zika), ii))
  
gf.list_plots(conditions, list_of_filtered_dfs, "Zika_454_29-10-19", 50)

#%% Ler PDBs

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

#%% Transformar sequencias em nós
#Dividir PDBs em listas de listas
pdbs_dict_list = fc.prepare_PDBs(pdbsNodesFiles)

#Dicinário de PDBs com listas de cadeias e nós
pdbs_dict_node_list = fc.dicts_of_pdbs(pdbs_dict_list)
#Lista de amostras em forma de nós
samples_node_list = fc.df_to_node_list(raw_df_aln)

get_ipython().magic("time")
#%% Alinhar
import Codes.smithWaterman as sw

sample_pdb = {"5GS6": pdbs_dict_node_list["5GS6"]}
sample_samples = samples_node_list[:10]

obj = sw.smithWaterman()
seqIds, seqPos, cut = obj.constructor(2, -1, -1, sample_samples[0], sample_pdb["5GS6"][0], True, False)



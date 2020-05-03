#%% Imports
import pandas as pd
from tqdm import tqdm

# My imports
import Codes.read as rd
import Codes.function as fc
import Codes.graph as gf
#%% Reading files

def read_files(path):
    #Leio o arquivo alinhado
    aux = path.split("/")
    path = aux[0]+"/"
    name = aux[-1]
    raw_aln_zika = rd.read_aligned_files(path, name)
    #Converto alinhado para um DF
    raw_df_aln = rd.seq_to_df(raw_aln_zika)
    
    return raw_df_aln

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

#%% Ploting


def plot_flow(path, conditions):
    
    #Criar
    raw_aln_zika = read_files(path)   
    counted_df = read_Or_create(raw_aln_zika)#, False)
    #Plotar
    list_of_filtered_dfs = []
    
    for ii in conditions:
      list_of_filtered_dfs.append(fc.filter_criteria(counted_df, len(raw_aln_zika), ii))
      
    gf.list_plots(conditions, list_of_filtered_dfs, "Zika_454_29-10-19", 50)

    return raw_aln_zika

conditions = [.5,1,2]
path = "read/"
name = "Zika Full 454 29-10-19.fasta"
GB_raw_df_aln = plot_flow(path+name, conditions)

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
pdbs_tuple_list = fc.prepare_PDBs(pdbsNodesFiles)
#Dicinário de PDBs com listas de cadeias e nós
pdbs_list_of_tuples = fc.pdbs_ids_to_nodes(pdbs_tuple_list)
#Lista de amostras em forma de nós
samples_node_list = fc.df_to_node_list(GB_raw_df_aln)

#%% Alinhar

import time
import Codes.smithWaterman as sw
start = time.time()

#Uma execução
# obj = sw.smithWaterman()
# seqIds, seqPos, cut = obj.constructor(2, -1, -1, samples_node_list[0], pdbs_list_of_tuples[0][1], False, False)
#Salvar resultados em uma tabela?

#jeito burro - Sequencial - 26.27 sec
# for sample in tqdm(sample_samples):
#     for key, pdb in sample_pdb.items():
#         for chain in pdb:
#             obj = sw.smithWaterman()
#             seqIds, seqPos, cut = obj.constructor(2, -1, -1, sample, chain, False, False)
     
    
end = time.time()
print(end - start, "sec")
#%%Salvar resultado em um pdb
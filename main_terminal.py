# -*- coding: utf-8 -*-
"""
Created on Sat May  2 16:39:11 2020

@author: ivana
"""
# Imports
import pandas as pd
from tqdm import tqdm
import sys
import argparse

# My imports
import Codes.read as rd
import Codes.function as fc
import Codes.graph as gf


def read_files(path):
	#path = "read/"
	#name = "Zika Full 454 29-10-19.fasta"
	aux = path.split("/")
	path = aux[0]+"/"
	name = aux[-1]
	raw_aln_zika = rd.read_aligned_files(path, name)
	#Converto alinhado para um DF
	raw_df_aln = rd.seq_to_df(raw_aln_zika)

	return raw_df_aln

def read_Or_create(raw_df_aln, read=True):
    if read:
        #lê salvo
        counted_df = pd.read_csv("Saved/1_Counting_df_zika_454_29-10-19X.csv")
    else:
        #Df vazio
        counted_df = pd.DataFrame()
        #função que faz a transposta de Seq e conta as ocorrencias
        fc.transpose_seq_and_count(raw_df_aln.Seq, counted_df)
        #Limpa os Nan e corrige a Pos
        counted_df = counted_df.fillna(0).reset_index().rename(columns={"index" : "Pos"})
        counted_df.Pos = counted_df.Pos.apply(lambda x: x+1)
        #Salva o trabalho
        counted_df.to_csv("Saved/1_Counting_df_zika_454_29-10-19X.csv", index=False)
    
    return counted_df

def flow():
	# function_name = argv.pop(0)
	# read_path = argv.pop(0)
	# filter = argv

	parser = argparse.ArgumentParser(description="Plot de SNPs")
	parser.add_argument("--filtro", "-f", required=True, nargs='+', type=float, help="Valores utilizados para o filtro")
	parser.add_argument("--path", "-p", required=True, help="Local do source")

	args = parser.parse_args()

	#Criar
	raw_df_aln = read_files(args.path)
	counted_df = read_Or_create(raw_df_aln, False)
	#Plotar 
	conditions = args.filtro
	list_of_filtered_dfs = []

	for ii in conditions:
		list_of_filtered_dfs.append(fc.filter_criteria(counted_df, len(raw_aln_zika), ii))
	
	gf.list_plots(conditions, list_of_filtered_dfs, "Zika_454_29-10-19X", 50)

	return 0

sys.exit(flow())
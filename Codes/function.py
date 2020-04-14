from tqdm import tqdm
import pandas as pd
import numpy as np
from collections import Counter 

def percentage_to_quantity(total, prc):
    """
    Converts percentage to a integer value
    Parameters
    ----------
    total : integer
        Total amount of samples
    prc : integer
        choosed percentage
    Returns
    -------
    integer
        Value equivalent of given percentage in the samples amount
    """
    return (prc*total)/100

def filter_criteria(df, total, condition=1):
    """
    Filtering data using pandas.
    Parameters
    ----------
    df : DataFrame
        DataFrame with the info
    total : integer
        Total of samples
    condition : integer, optional
        Criteria for the filter. The default is 1.
    Returns
    -------
    temp : DataFrame
        Filtered data
    """
    temp = df[df < total - percentage_to_quantity(total, condition)]
    temp.Pos = df.Pos
      
    return temp

def transpose_seq_and_count(column, counted_df):
    """
    Faz a transposta da coluna Seq e contabiliza as repetições por posição
    nas sequências

    Parameters
    ----------
    column : Pandas.Series
        A coluna Seq, com todas as amostras que serão contados
    counted_df : pd.DataFrame
        Dtaframe vazio onde será armazenado o resultado

    Returns
    -------
    None.

    """
    #Deixa todas as listas em Upper Case
    aux = column.apply(str.upper)
    #Transformei cada seq.row em lista
    aux = aux.apply(list)
    #Tranformei a coluna seq em uma lista de listas e fiz a Transposta
    aux = np.transpose(list(aux))
    
    #Paraca cada posição(lista) contar Counter
    for idx, row in tqdm(enumerate(aux), total=len(aux)):
        row = Counter(row)
        for nuc, count in row.items():
            if nuc == "U":
                nuc = "T"
            counted_df.loc[idx, nuc] = count
            
def prepare_PDBs(pdbs_dict):
    """
    Recebe dicionário com pdbs no formato de DataFrame e retorna um dicionário
    com uma lista de listas, cada lista contém od pdbs divididos em cadeia

    Parameters
    ----------
    pdbs_dict : Dict
        Dicionário de DataFrames usado na entrada

    Returns
    -------
    result : Dict
        dicionário de listas de listas

    """
    result = {}
    for key, value in pdbs_dict.items():
        #divide o df por cadeias
        temp_chain_list = value.Chain.unique().tolist()
        aux = []
        for chain in temp_chain_list:
            #filtra cadeias e salva em lista
            filtered_df = value[value.Chain == chain]
            aux.append(list(zip(filtered_df.NodeId, filtered_df.Degree)))
        result[key] = aux
        
    return result
from tqdm import tqdm
import numpy as np
from collections import Counter
from Bio.SeqUtils import seq1
#MyImports
import Codes.node as nd
import Codes.smithWaterman as sw

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
    Recebe dicionário com pdbs no formato de DataFrame e retorna uma lista
    de tuplas, com o nome do PDB e uma lista com os os ids dos nós 

    Parameters
    ----------
    pdbs_dict : Dict
        Dicionário de DataFrames usado na entrada

    Returns
    -------
    result : list
        lista de tuplas -> [(name, [node1, node2,...]), (name, [node1, node2,...])]

    """
    result = []
    for key, df in pdbs_dict.items():
        #divide o df por cadeias
        chains_on_df = df.Chain.unique().tolist()
        for chain in chains_on_df:
            #filtra cadeias e salva em lista
            filtered_df = df[df.Chain == chain]
            result.append((key, list(zip(filtered_df.NodeId, filtered_df.Degree))))
        
    return result

def nodeIds_to_node(id):
    """
     transformar elmentos das listas em nós   

    Parameters
    ----------
    id : tuple
        tupla contendo o ID e o grau da amostra

    Returns
    -------
    Node
        Classe criada por mim, objeto Node -> ver node.py

    """
    splited_Id = id[0].split(":")
    degree = id[1]
    return nd.node(seq1(splited_Id[3]), degree, id[0], splited_Id[1], 0)

def list_of_nodeIds(tuple_list):
    """
    Chama map para elementos da sublista

    Parameters
    ----------
    tuple_list : list
        lista de tuplas, NodeId e Degree

    Returns
    -------
    sub_list : list
        lista de Nodes -> nodes.py

    """
    return list(map(nodeIds_to_node, tuple_list[1]))

def pdbs_ids_to_nodes(tuple_list):
    """
    Chama map para lista de listas

    Parameters
    ----------
    tuple_list : list
        listas de tuplas
        0: PDB name
        1: lista de NodeIds

    Returns
    -------
    pdbs_dict : list
        lista de tuplas -> [(NodeId, [nodes list]),(NodeId,...]

    """
    return [(tp[0], list_of_nodeIds(tp)) for tp in tuple_list]


def sample_tuple_to_nodes(tuple):
    return nd.node(tuple[0], 0, tuple[2], 0, tuple[1])
  
def sample_to_tuple_to_nodes(seq):
    name = seq[1]
    seq = seq[0]
    """
    Transforma Seq em uma lista de tuplas (Seq, Pos)

    Parameters
    ----------
    seq : String
        Sequencia de DNA da amostra

    Returns
    -------
    list
        Lita de tuplas, Nodes -> node.py

    """
    size = len(seq)
    temp = list(zip(seq, list(range(size)), [name]*size))
    return list(map(sample_tuple_to_nodes, temp))


def df_to_node_list(df):
    """
    Transformar df de samples em listas de tuplas

    Parameters
    ----------
    df : DataFrame
        DataFRame com amostras

    Returns
    -------
    list
        Listas de tuplas.

    """
    return list(map(sample_to_tuple_to_nodes, list(zip(df.Seq, df.Name))))

def one_step(sample, pdb_chain):
    """
    Alinha duas sequencias.

    Parameters
    ----------
    sample : list
        lista de nós da sequencia da amostra
    pdb_chain : list
        lista de nós da sequencia do PDB

    Returns
    -------
    seqIds : list
        Lista dos ids do pdb
    seqPos : list
        lista das posições da amosta que alinhou contra o pdb
    cut : int
        Volor de cobertura o alinhamento

    """
    #print(len(pdb_chain), type(pdb_chain), pdb_chain)
    obj = sw.smithWaterman()
    seqIds, seqPos, cover = obj.constructor(2, -1, -1, sample, pdb_chain, True, False)
    
    return ",".join(seqIds), ",".join(seqPos), cover
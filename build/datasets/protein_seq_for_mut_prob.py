"""
Generate fasta proteins sequences to be backtranslated to DNA and initialize a dataframe
that will be used to store Gene/Uniprot ID with protein and DNA sequences.

Given a directory including PDB structures (tested on structures from AF DB), generate a series 
of fasta files including the protein sequences. These files will be used in the next step to 
retrieve DNA sequence by backtranslation. It also generate a dataframe that include the protein 
sequence of each protein included in the input file. This dataframe will be annotated with DNA 
sequences in the next step.

###################################### EXAMPLE USAGE ####################################

python3 protein_seq_for_mut_prob.py -i /workspace/datasets/alphafold_features/AF_homo_sapiens_pred/ \
-o ../../datasets/temp/fasta_seq -O ../../datasets/seq_for_mut_prob.csv \
-u ../../datasets/af_uniprot_to_gene_id.json
      
#########################################################################################
"""

import os
import json
import pandas as pd
import json
import argparse
import warnings
from progressbar import progressbar
from utils import uniprot_to_hugo, get_pdb_path_list_from_dir, get_af_id_from_pdb, get_seq_from_pdb


def get_seq_df_from_dir(input_path, uniprot_to_gene_dict):
    """
    Parse any PDB structure from a given directory and create 
    a dataframe including HUGO symbol, Uniprot-ID, AF fragment,
    and protein sequence.
    """
    
    # Get all PDB path in directory
    list_prot_path = get_pdb_path_list_from_dir(input_path)

    # Get Uniprot ID, HUGO, F and protein sequence of any PDB in dir
    gene_lst = []
    uni_id_lst = []
    f_lst = []
    seq_lst = []
    id_not_found = []

    for path_structure in progressbar(list_prot_path):
        identifier = get_af_id_from_pdb(path_structure)
        uniprot_id, f = identifier.split("-F")

        if uniprot_id in uniprot_to_gene_dict.keys():
            gene = uniprot_to_gene_dict[uniprot_id]
            seq = "".join(list(get_seq_from_pdb(path_structure)))
            gene_lst.append(gene)
            uni_id_lst.append(uniprot_id)
            f_lst.append(f)
            seq_lst.append(seq)
        else:
            id_not_found.append(uniprot_id)                         

    if len(id_not_found) > 0:                                  ## I could output a df with IDs not found
        warnings.warn(f"WARNING........... {len(id_not_found)} Uniprot-ID not found in the Uniprot-HUGO mapping dictionary")

    seq_df = pd.DataFrame({"Gene" : gene_lst, "Uniprot_ID" : uni_id_lst, "F" : f_lst, "Seq" : seq_lst}).sort_values(["Gene", "F"])

    return seq_df


def protein_df_to_fasta(seq_df, output_dir, max_seq=500, max_size=900000, v=True):
    """
    From a dataframe including genes, proteins IDs, and 
    amino acid sequences, generate one or more fasta files 
    including protein sequences to be backtranslated to DNA
    by https://www.ebi.ac.uk/Tools/st/emboss_backtranseq/
    
    Arguments
    ---------
    seq_df: pandas df
        It must include the following cols: "Gene", "Uniprot_ID", "Seq"
    max_seq: int
        Max number of sequences in a single fasta file
    max_size: int
        Max bytes of a single fasta file
    fasta_out: str
    
    Returns
    -------
    seq_df: pandas df
        The dataframe include two additional columns: "Batch" indicates
        the fasta file where the protein seq was saved and "N_prot" is 
        the protein number in the file.
    """
    
    i = 1
    n = 1
    filename = f"proteins/af_proteins{n}.fasta"
    file = open(f"{output_dir}/{filename}", 'w')

    # Iterate through the AF proteins seq
    for j, row in progressbar(seq_df.iterrows()):
        
        # Write tlo file
        size = os.path.getsize(f"{output_dir}/{filename}") 
        if i > max_seq or size > max_size:
            if v:
                print(f"File {output_dir}{filename}")
                print(f"Batch {n}, {size} bytes, {i-1} sequences\n")
            file.close()    
            i = 1                                                # Reset protein # in batch 
            n += 1                                               # Add 1 to batch counter
           
            filename = f"proteins/af_proteins{n}.fasta"
            file = open(f"{output_dir}/{filename}", 'w')
            
        gene = row["Gene"]
        prot_id = row["Uniprot_ID"]
        seq = row["Seq"]
        seq_df.at[j, "Batch"] = n
        seq_df.at[j, "N_prot"] = i    
        
        if pd.isnull(gene):
            gene = "NA"
        if pd.isnull(prot_id):
            gene = "NA"
        file.write(f">{gene}|{prot_id}|{format(i, '03d')}\n{seq}\n")
        
        i += 1                                                   # Add 1 to protein # in batch 
        
    if v:
        print(f"File {output_dir}{filename}")
        print(f"Batch {n}, {size} bytes, {i-1} sequences\n")
    file.close()
    
    # Write protein sequence file to retrieve all info 
    seq_df.Batch = pd.to_numeric(seq_df.Batch, downcast = "integer")
    seq_df.N_prot = pd.to_numeric(seq_df.N_prot, downcast = "integer")
    
    return seq_df


def main():
    
    ## Parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Path to directory including PDB structures", type=str, required=True)
    parser.add_argument("-o", "--output_dir_fasta", help="Output directory to save fasta sequences", type=str, 
                        default="../../datasets/temp/fasta_seq/")     
    parser.add_argument("-O", "--output_seq_df", help="Output path to save the dataframe including all sequences", type=str, 
                        default="../../datasets/seq_for_mut_prob.csv")                                                   
    parser.add_argument("-v", "--verbose", help="Verbose", type=int, default=0)          
    parser.add_argument("-u", "--uniprot_to_gene_dict", help="Path to a dictionary including Uniprot_ID : HUGO symbol mapping", type=str, 
                        default="../../datasets/af_uniprot_to_gene_id.json")  
    parser.add_argument("-j", "--join_fragments", help="if 1 join fragments of each protein, if 0 return them as unique elements", type=int, default=0)  

    args = parser.parse_args()
    input = args.input
    uniprot_to_gene_dict = args.uniprot_to_gene_dict
    output_dir_fasta = args.output_dir_fasta
    output_seq_df = args.output_seq_df
    join_fragments = args.join_fragments
    verbose = args.verbose


    ## Initialize

    print(f"Join fragments of large proteins: {join_fragments}")

    # Create necessary folders
    proteins_path = f"{output_dir_fasta}/proteins"
    dna_path = f"{output_dir_fasta}/dna"
    if not os.path.exists(proteins_path):
        os.makedirs(proteins_path)
    if not os.path.exists(dna_path): 
        os.makedirs(dna_path)

    # Load Uniprot ID to HUGO mapping
    # neighbours_df = pd.read_pickle(neighbours_path)
    if uniprot_to_gene_dict is not None:
        uniprot_to_gene_dict = json.load(open(uniprot_to_gene_dict)) 
    else:
        uniprot_to_gene_dict = uniprot_to_hugo()  

    # Create a dataframe with protein sequences
    print("\nGenerating sequence df..")
    seq_df = get_seq_df_from_dir(input, uniprot_to_gene_dict)

    # Create fasta files with protein sequences
    print(f"\nSaving fasta files for {len(seq_df)} proteins..")
    seq_df = protein_df_to_fasta(seq_df, output_dir=output_dir_fasta, v=verbose)
    seq_df.to_csv(output_seq_df, index=False)
    print(f"Dataframe including sequences is saved in: {output_seq_df}")
    print(f"Fasta files are saved in: {output_dir_fasta}")


if __name__ == "__main__":
    main()
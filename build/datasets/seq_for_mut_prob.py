"""
Generate a pandas dataframe including identifiers mapped to protein and DNA sequences.

Given a directory including PDB structures (tested on structures from AF DB), extract all
protein sequences, use EMBOSS backtranseq back translate proteins sequences into DNA, and
generate a dataframe including HUGO symbol, Uniprot_ID, protein, and DNA sequences.
This dataframe is required to get the probability of each residue to mutate based on the 
mutation profile (mutation rate in trinucleotide contexts) of the cohort. Which is then used
to get the probability of a certain volume to be hit by a missense mutation.

###################################### EXAMPLE USAGE ####################################

python3 seq_for_mut_prob.py -i ../../datasets/pdb_structures/ \
-o ../../datasets/seq_for_mut_prob.csv \
-u ../../datasets/af_uniprot_to_gene_id.json \ 
-s "Mus musculus"
      
#########################################################################################
"""

import json
import pandas as pd
import json
import argparse
import warnings
from progressbar import progressbar
import requests
import time
import os
import re
from utils import uniprot_to_hugo, get_pdb_path_list_from_dir, get_af_id_from_pdb, get_seq_from_pdb, get_seq_similarity, translate_dna


def get_seq_df_from_dir(input_path, uniprot_to_gene_dict):
    """
    Parse any PDB structure from a given directory and create 
    a dataframe including HUGO symbol, Uniprot-ID, AF fragment,
    and protein sequence.
    """
    
    # Get all PDB path in directory whose IDs can be mapped to HUGO symbols
    list_prot_path = get_pdb_path_list_from_dir(input_path)
    list_prot_path = [path for path in list_prot_path if get_af_id_from_pdb(path).split("-F")[0] in uniprot_to_gene_dict.keys()]
    pdb_not_mapped = set([get_af_id_from_pdb(path).split("-F")[0] for path in list_prot_path if get_af_id_from_pdb(path).split("-F")[0] not in uniprot_to_gene_dict.keys()])
    if len(pdb_not_mapped) > 0:                                 
        warnings.warn(f"{len(pdb_not_mapped)} Uniprot-ID not found in the Uniprot-HUGO mapping dictionary")

    # Get Uniprot ID, HUGO, F and protein sequence of any PDB in dir
    gene_lst = []
    uni_id_lst = []
    f_lst = []
    seq_lst = []

    for path_structure in progressbar(list_prot_path):
        uniprot_id, f = get_af_id_from_pdb(path_structure).split("-F")
        gene = uniprot_to_gene_dict[uniprot_id]
        seq = "".join(list(get_seq_from_pdb(path_structure)))
        gene_lst.append(gene)
        uni_id_lst.append(uniprot_id)
        f_lst.append(f)
        seq_lst.append(seq)

    seq_df = pd.DataFrame({"Gene" : gene_lst, "Uniprot_ID" : uni_id_lst, "F" : f_lst, "Seq" : seq_lst}).sort_values(["Gene", "F"])

    return seq_df


def backtranseq(protein_seqs, organism = "Homo sapiens"):
    """
    Perform backtranslation from proteins to 
    DNA sequences using EMBOS backtranseq.    
    """
    
    # Define the API endpoints
    run_url = "https://www.ebi.ac.uk/Tools/services/rest/emboss_backtranseq/run"
    status_url = "https://www.ebi.ac.uk/Tools/services/rest/emboss_backtranseq/status/"
    result_url = "https://www.ebi.ac.uk/Tools/services/rest/emboss_backtranseq/result/"
    
    # Define the parameters for the API request (an email address must be included)
    params = {"email": "emaple.email@irbbarcelona.org",
              "sequence": protein_seqs,
              "outseqformat": "plain",
              "molecule": "dna",
              "organism": organism}
    
    # Submit the job request and retrieve the job ID
    response = requests.post(run_url, data=params)
    job_id = response.text.strip()

    # Wait for the job to complete
    status = ""
    while status != "FINISHED":
        time.sleep(10)
        result = requests.get(status_url + job_id)
        status = result.text.strip()

    # Retrieve the results of the job
    result = requests.get(result_url + job_id + "/out")
    dna_seq = result.text.strip()
    
    return dna_seq


def batch_backtranseq(df, batch_size, organism = "Homo sapiens"):
    """
    Given a dataframe including protein sequences, it divides the 
    sequences into batches of a given size and run EMBOSS backtranseq 
    (https://www.ebi.ac.uk/Tools/st/emboss_backtranseq/) to translate 
    them into DNA sequences.
    """
    
    batches = df.groupby(df.index // batch_size)
    lst_batches = []
        
    # Iterate over batches
    for _, batch in progressbar(batches):

        # Get input format for backtranseq 
        batch_seq = "\n".join(batch.reset_index(drop=True).apply(lambda x: f'>\n{x["Seq"]}', axis=1).values)

        # Run backtranseq
        batch_dna = backtranseq(batch_seq, organism = organism)

        # Parse output
        batch_dna = re.split(">EMBOSS_\d+", batch_dna.replace("\n", ""))[1:]

        batch["Seq_dna"] = batch_dna
        lst_batches.append(batch)
        
    return pd.concat(lst_batches)


def add_extra_genes_to_seq_df(seq_df, uniprot_to_gene_dict):
    """
    If multiple genes are mapping to a given Uniprot_ID, add 
    each gene name with corresponding sequence info to the seq_df.
    """
    
    lst_extra_genes = []

    for _, seq_row in seq_df.iterrows():

        uni_id = seq_row["Uniprot_ID"]
        gene_id = uniprot_to_gene_dict[uni_id]
        if pd.isnull(gene_id) == False:

            gene_id = gene_id.split(" ")
            if len(gene_id) > 1:

                for gene in gene_id:
                    if gene != seq_row["Gene"]:
                        lst_extra_genes.append((gene, seq_row["Uniprot_ID"], seq_row["F"], seq_row["Seq"], seq_row["Seq_dna"]))

    seq_df_extra_genes = pd.DataFrame(lst_extra_genes, columns=["Gene", "Uniprot_ID", "F", "Seq", "Seq_dna"])
    
    return pd.concat((seq_df, seq_df_extra_genes))


def main():
    
    ## Parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Path to directory including PDB structures", type=str, required=True)
    parser.add_argument("-o", "--output_seq_df", help="Output path to save the dataframe including all sequences", type=str, 
                        default="../../datasets/seq_for_mut_prob.csv")                                                       
    parser.add_argument("-u", "--uniprot_to_gene_dict", help="Path to a dictionary including Uniprot_ID : HUGO symbol mapping", type=str)  
    parser.add_argument("-s", "--organism", help="Binominal nomenclature of the organism", type=str, default="Homo sapiens") 

    args = parser.parse_args()
    input = args.input
    uniprot_to_gene_dict = args.uniprot_to_gene_dict
    output_seq_df = args.output_seq_df
    organism = args.organism

    # Load Uniprot ID to HUGO mapping
    if uniprot_to_gene_dict is not None:
        uniprot_to_gene_dict = json.load(open(uniprot_to_gene_dict)) 
    else:
        print("Retrieving Uniprot_ID to Hugo symbol mapping information..")
        uniprot_ids = os.listdir(input)
        uniprot_ids = [uni_id.split("-")[1] for uni_id in list(set(uniprot_ids)) if uni_id.endswith(".pdb")]
        uniprot_to_gene_dict = uniprot_to_hugo(uniprot_ids)  

    # Create a dataframe with protein sequences
    print("\nGenerating sequence df..")
    seq_df = get_seq_df_from_dir(input, uniprot_to_gene_dict)

    # Annotate df with DNA sequences
    print("\nPerforming back translation")
    seq_df = batch_backtranseq(seq_df, 500, organism=organism)
    
    # Add multiple genes mapping to the same Uniprot_ID
    seq_df = add_extra_genes_to_seq_df(seq_df, uniprot_to_gene_dict)
    seq_df = seq_df.dropna(subset=["Gene"])
    
    # Save and assess similarity
    seq_df.to_csv(output_seq_df, index=False)
    sim_ratio = sum(seq_df.apply(lambda x: get_seq_similarity(x.Seq, translate_dna(x.Seq_dna)), axis=1)) / len(seq_df)
    if sim_ratio < 1:                       
        warnings.warn(f"WARNING........... some problems occurred during back translation, protein and translated DNA similarity: {sim_ratio}")
    print(f"Dataframe including sequences is saved in: {output_seq_df}")


if __name__ == "__main__":
    main()
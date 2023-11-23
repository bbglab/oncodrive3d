"""
Contain functions to compute the per-residues probability of missense 
mutation of any protein given the mutation profile of the cohort.
"""


import argparse
import json
import os
from itertools import product

import daiquiri
import numpy as np
import pandas as pd
from progressbar import progressbar

from mutability import Mutabilities
from mutability import init_mutabilities_module

logger = daiquiri.getLogger(__logger_name__ + ".utils.miss_mut_prob")


def get_unif_gene_miss_prob(size):
    """
    Get a uniformly distributed gene missense mutation 
    probability vector.
    """
    
    vector = np.ones(size)
    vector[0] = 0
    vector[-1] = 0
    
    return vector / sum(vector)


def mut_rate_vec_to_dict(mut_rate, v=False):
    """
    Convert the vector of mut mut_rate of 96 channels to a dictionary of 192 
    items: the keys are mutations in trinucleotide context (e.g., "ACA>A") 
    and values are the corresponding mut rate (frequency of mut normalized 
    for the nucleotide content).
    """
    
    cb  = dict(zip('ACGT', 'TGCA'))
    if v: logger.debug("Mut1\tMut2\tN_context\tMut_rate")
    mut_rate_dict = {}
    i = 0
    for ref in ['C', 'T']:
        for alt in cb.keys():
            if ref == alt:
                continue
            else:
                for p in product(cb.keys(), repeat=2):
                    mut = f"{p[0]}{ref}{p[1]}>{alt}"
                    cmut = f"{cb[p[1]]}{cb[ref]}{cb[p[0]]}>{cb[alt]}"
                    mut_rate_dict[mut] = mut_rate[i]
                    mut_rate_dict[cmut] = mut_rate[i]
                    if v: logger.debug(f"{mut}\t{cmut}\t{i}\t\t{mut_rate[i]}")
                    i +=1
                    
    return mut_rate_dict


def get_codons(dna_seq):
    """
    Get the list of codons from a DNA sequence.
    """
    
    return [dna_seq[i:i+3] for i in [n*3 for n in range(int(len(dna_seq) / 3))]]


def translate_dna_to_prot(dna_seq, gencode):
    """
    Translate a DNA sequence into amino acid sequence.
    """
    
    return "".join([gencode[codon] for codon in get_codons(dna_seq)])


def colored(text, r=255, g=0, b=0):
    
    start = f"\033[38;2;{r};{g};{b}m"
    end = "\033[38;2;255;255;255;0m"
    
    return f"{start}{text}{end}"


def get_miss_mut_prob(dna_seq, mut_rate_dict, mutability=False, get_probability=True, v=False):
    """
    Generate a list including the probabilities that the 
    codons can mutate resulting into a missense mutations.
    
    Arguments
    ---------
    dna_seq: str
        Sequence of DNA
    mut_rate_dict: dict
        Mutation rate probability as values and the 96 possible
        trinucleotide contexts as keys
    gencode: dict
        Nucleotide as values and codons as keys
        
    Returns
    -------
    missense_prob_vec: list
        List of probabilities (one for each codon or prot res) 
        of a missense mutation  
    """

    # Initialize
    gencode = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}

    # Get all codons of the seq
    codons = get_codons(dna_seq)

    # Iterate through all codons except the first and last
    missense_prob_vec = [0]                                 # Assign 0 prob to the first residue
    
    # TODO revise whether these needs to be len(codons)-1 or len(codons) since the range already does -1
    # we discussed on 2023-11-21 and this -1 in the lenght of the codons could be removed
    # but we should check if the trinucl2 would still work, it will probably raise an error because
    # it will not find the nucleotide corresponding to the first position of the next codon
    for c in range(1, len(codons)-1):
        missense_prob = 0
        codon = codons[c]
        trinucl0 = f"{codons[c-1][2]}{codons[c][0:2]}"
        trinucl1 = f"{codons[c]}"
        trinucl2 = f"{codons[c][1:3]}{codons[c+1][0]}"

        # Print codon info
        if v:
            logger.debug(f"\n{colored('>>>')} Codon n: {c}", "\tRef AA:", aa, "\tCodon:", colored(codon))
            logger.debug("")
            logger.debug(f"\t\t{codons[c-1]}{colored(codons[c])}{codons[c+1]}")
            logger.debug(f"\t\t..{colored(trinucl0, 0, 0, 255)}....")
            logger.debug(f"\t\t...{colored(trinucl1, 0, 0, 255)}...")
            logger.debug(f"\t\t....{colored(trinucl2, 0, 0, 255)}..")
            logger.debug("")

        # Iterate through the possible contexts of a missense mut
        for i, trinucl in enumerate([trinucl0, trinucl1, trinucl2]):
            ref = trinucl[1]
            aa = gencode[codon]
            if v: 
                logger.debug(">> Context:", colored(trinucl, 0, 0, 255), "\n   Ref:", ref, "\n")
                logger.debug("Mut      Mut_prob                Alt_codon   Alt_AA")

            # Iterate through the possible alt 
            for alt in "ACGT":
                if alt != ref:         
                    alt_codon = [n for n in codon]
                    alt_codon[i] = alt
                    alt_codon = "".join(alt_codon)
                    alt_aa = gencode[alt_codon]  
                    # If there is a missense mut, get prob from context and sum it
                    if alt_aa != aa and alt_aa != "_":
                        if not mutability:
                            mut = f"{trinucl}>{alt}"    # query using only the trinucleotide change
                            if v: print(f"{trinucl}>{alt}", "\t", mut_rate_dict[mut], "\t", alt_codon, "\t    ", alt_aa, )
                            if mut in mut_rate_dict:
                                missense_prob += mut_rate_dict[mut]
                            else:
                                missense_prob += 0

                        else:
                            # TODO this has not been tested
                            cdna_pos = (c * 3) + i  # compute the cDNA position of the residue
                            if cdna_pos in mut_rate_dict:
                                missense_prob += mut_rate_dict[cdna_pos].get(alt, 0)
                            else:
                                missense_prob += 0


            if v: logger.debug("")

        if v: logger.debug(f">> Prob of missense mut: {missense_prob:.3}\n")
        missense_prob_vec.append(missense_prob)
    missense_prob_vec.append(0)                               # Assign 0 prob to the last residue

    # Convert into probabilities
    if get_probability:
        missense_prob_vec = np.array(missense_prob_vec) / sum(missense_prob_vec)
    
    return list(missense_prob_vec)



def get_miss_mut_prob_dict(mut_rate_dict, seq_df, mutability=False, mutability_config=None, v=False):
    """
    Given a dictionary of mut rate in 96 contexts (mut profile) and a 
    dataframe including Uniprot ID, HUGO symbol and DNA sequences, 
    get a dictionary with UniprotID-Fragment as keys and corresponding 
    vectors of missense mutation probabilities as values.
    """

    miss_prob_dict = {}

    if mutability:
        # TODO if the execution time of this step is too long we could
        # parallelize all these loops so that each gene is done in parallel

        # Process any Protein/fragment in the sequence df
        if "F" in seq_df.columns:
            for _, row in seq_df.iterrows():
                # Mutabilities
                mutability_dict = Mutabilities(row.Uniprot_ID, row.Chr, row.Exons_coord, len(row.Seq_dna), row.reverse_strand, mutability_config).mutabilities_by_pos
                miss_prob_dict[f"{row.Uniprot_ID}-F{row.F}"] = get_miss_mut_prob(row.Seq_dna, mutability_dict, mutability=True, v=v)
        else:
            for _, row in seq_df.iterrows():
                mutability_dict = Mutabilities(row.Uniprot_ID, row.Chr, row.Exons_coord, len(row.Seq_dna), row.reverse_strand, mutability_config).mutabilities_by_pos
                miss_prob_dict[f"{row.Uniprot_ID}"] = get_miss_mut_prob(row.Seq_dna, mutability_dict, mutability=True, v=v)

    else:
        # Process any Protein/fragment in the sequence df
        if "F" in seq_df.columns:
            for _, row in seq_df.iterrows():
                miss_prob_dict[f"{row.Uniprot_ID}-F{row.F}"] = get_miss_mut_prob(row.Seq_dna, mut_rate_dict, v=v)
        else:
            for _, row in seq_df.iterrows():
                miss_prob_dict[f"{row.Uniprot_ID}"] = get_miss_mut_prob(row.Seq_dna, mut_rate_dict, v=v)

    return miss_prob_dict


# TODO revise if this main() function should be here or not
def main():
    
    ## Parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", 
                        help="Input path to the directory including the mutation profile of the cohorts (list of 96 floats or dict of 192 items)", 
                        type=str, required=True)
    parser.add_argument("-o", "--output", 
                        help="Output path to save the dictionaries of missense mutation probability of each protein of the cohorts", 
                        type=str, required=True)                        
    parser.add_argument("-s", "--seq_df", 
                        help="Path to the dataframe including DNA and protein seq of all gene/proteins (all AF predicted ones)", 
                        type=str, 
                        default="../../datasets/seq_for_mut_prob.csv")                                         
    parser.add_argument("-v", "--verbose", help="Verbose", type=int, default=0)          

    args = parser.parse_args()
    input_path = args.input
    output_path = args.output
    path_seq_df = args.seq_df
    verbose = args.verbose


    logger.info("Input directory:", input_path)
    logger.info("Output directory:", output_path)
    logger.info("Path to DNA sequences df", path_seq_df, "\n")

    # Create necessary folder
    if not os.path.exists(output_path):
        os.makedirs(output_path)


    # Iterate through all mut profiles path in the directory
    path_profiles = [f"{input_path}/{file}" for file in os.listdir(input_path) if file.endswith(".json")]
    for path_mut_profile in progressbar(path_profiles):

        # Filename
        filename = path_mut_profile.split("/")
        filename = filename[len(filename)-1]
        filename = filename.split(".json")[0]

        # Load mut profile (mut rate) and convert into dictionary
        with open(path_mut_profile) as json_file:
            mut_rate_dict = json.load(json_file)
        if not isinstance(mut_rate_dict, dict):
            mut_rate_dict = mut_rate_vec_to_dict(mut_rate_dict)

        # Get the per-residue miss mut prob for each protein and add it to a dict
        seq_df = pd.read_csv(path_seq_df)
        miss_prob_dict = {}
        logger.info(f" Processing {len(seq_df)} proteins/fragments in {filename}..")
        miss_prob_dict = get_miss_mut_prob_dict(mut_rate_dict=mut_rate_dict, seq_df=seq_df, v=verbose)

        # Save
        with open(f"{output_path}/{filename}.miss_mut_prob.json", "w") as fp:
            json.dump(miss_prob_dict, fp)


if __name__ == "__main__":
    # main()
    mutab_config = json.load(open('/home/fcalvet/Documents/dev/clustering_3d/test/normal_tests/mutability_config.json'))
    init_mutabilities_module(mutab_config)

    seq_df = pd.read_csv("/home/fcalvet/Documents/dev/clustering_3d/test/normal_tests/seq_df_toy_coord.csv", header = 0)
    seq_df["Exons_coord"] = [ eval(exon_coord) for exon_coord in seq_df['Exons_coord'] ]

    miss_prob_dict = get_miss_mut_prob_dict(mut_rate_dict=None, seq_df=seq_df,
                                            mutability=True, mutability_config=mutab_config)
    print(miss_prob_dict)
    
#     chrom = 17
#     # exons = eval("[(7676594, 7676521)]")
#     # seq_len = len("ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACT")
#     exons = eval("[(7676594, 7676521), (7676403, 7676382), (7676272, 7675994), (7675236, 7675053), (7674971, 7674859), (7674290, 7674181), (7673837, 7673701), (7673608, 7673535), (7670715, 7670609), (7669690, 7669612)]")
#     seq_len = len("ATGGAGGAGCCCCAGAGCGACCCCAGCGTGGAGCCCCCCCTGAGCCAGGAGACCTTCAGCGACCTGTGGAAGCTGCTGCCCGAGAACAACGTGCTGAGCCCCCTGCCCAGCCAGGCCATGGACGACCTGATGCTGAGCCCCGACGACATCGAGCAGTGGTTCACCGAGGACCCCGGCCCCGACGAGGCCCCCAGGATGCCCGAGGCCGCCCCCCCCGTGGCCCCCGCCCCCGCCGCCCCCACCCCCGCCGCCCCCGCCCCCGCCCCCAGCTGGCCCCTGAGCAGCAGCGTGCCCAGCCAGAAGACCTACCAGGGCAGCTACGGCTTCAGGCTGGGCTTCCTGCACAGCGGCACCGCCAAGAGCGTGACCTGCACCTACAGCCCCGCCCTGAACAAGATGTTCTGCCAGCTGGCCAAGACCTGCCCCGTGCAGCTGTGGGTGGACAGCACCCCCCCCCCCGGCACCAGGGTGAGGGCCATGGCCATCTACAAGCAGAGCCAGCACATGACCGAGGTGGTGAGGAGGTGCCCCCACCACGAGAGGTGCAGCGACAGCGACGGCCTGGCCCCCCCCCAGCACCTGATCAGGGTGGAGGGCAACCTGAGGGTGGAGTACCTGGACGACAGGAACACCTTCAGGCACAGCGTGGTGGTGCCCTACGAGCCCCCCGAGGTGGGCAGCGACTGCACCACCATCCACTACAACTACATGTGCAACAGCAGCTGCATGGGCGGCATGAACAGGAGGCCCATCCTGACCATCATCACCCTGGAGGACAGCAGCGGCAACCTGCTGGGCAGGAACAGCTTCGAGGTGAGGGTGTGCGCCTGCCCCGGCAGGGACAGGAGGACCGAGGAGGAGAACCTGAGGAAGAAGGGCGAGCCCCACCACGAGCTGCCCCCCGGCAGCACCAAGAGGGCCCTGCCCAACAACACCAGCAGCAGCCCCCAGCCCAAGAAGAAGCCCCTGGACGGCGAGTACTTCACCCTGCAGATCAGGGGCAGGGAGAGGTTCGAGATGTTCAGGGAGCTGAACGAGGCCCTGGAGCTGAAGGACGCCCAGGCCGGCAAGGAGCCCGGCGGCAGCAGGGCCCACAGCAGCCACCTGAAGAGCAAGAAGGGCCAGAGCACCAGCAGGCACAAGAAGCTGATGTTCAAGACCGAGGGCCCCGACAGCGAC")

#     tot_s_ex = 0
#     for s, e in exons:
#         tot_s_ex += np.sqrt((e-s)**2) + 1
#     #print(tot_s_ex)

#     mutability_obj = Mutabilities("TP53", chrom, exons, seq_len, True, mutab_config)
    
#     # for s, e in exons:
#     #     for ii in range(min(s, e), max(s, e)+1):
#     #         print(ii)

# #    for key in mutability_obj.mutabilities_by_pos.keys():
# #        print(key)
# #        if len(mutability_obj.mutabilities_by_pos[key]) != 3:
# #            print(mutability_obj.mutabilities_by_pos[key])

#     for key in sorted(mutability_obj.mutabilities_by_pos):
#         # print(key)
#         print(key, mutability_obj.mutabilities_by_pos[key])
#         # if len(mutability_obj.mutabilities_by_pos[key]) != 3:
#         #     print(mutability_obj.mutabilities_by_pos[key])

#     # print(len(mutability_obj.mutabilities_by_pos))
#     # print(seq_len)
#     mutability_dict = mutability_obj.mutabilities_by_pos
    
#     # TODO raise an error here, well not here but within the load mutabilities function maybe?
#     if len(mutability_dict) != seq_len:
#         print("error")
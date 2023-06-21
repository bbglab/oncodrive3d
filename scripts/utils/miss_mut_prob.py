
"""
Process all mutation profiles in a given directory and, for each one, 
generate a dictionary having as keys UniprotID-Fragment (eg., P04217-F1)
and as values the corresponding per-residue miss mutation probability.

###################################### EXAMPLE USAGE ############################################

python3 miss_mut_prob.py -i ../../required_files/extra/mut_profile/ \
-o ../../required_files/extra/missense_mut_prob/ \
-s ../../required_files/seq_for_mut_prob.csv 

#################################################################################################

"""


import json
import pandas as pd
import numpy as np
import json
import argparse
from itertools import product
import os
from progressbar import progressbar


def mut_rate_vec_to_dict(mut_rate, v=False):
    """
    Convert the vector of mut mut_rate of 96 channels to a dictionary of 192 
    items: the keys are mutations in trinucleotide context (e.g., "ACA>A") 
    and values are the corresponding mut rate (frequency of mut normalized 
    for the nucleotide content).
    """
    
    cb  = dict(zip('ACGT', 'TGCA'))
    if v: print("Mut1\tMut2\tN_context\tMut_rate")
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
                    if v: print(f"{mut}\t{cmut}\t{i}\t\t{mut_rate[i]}")
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


def get_miss_mut_prob(dna_seq, mut_rate_dict, get_probability=True, v=False):
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
    for c in range(1, len(codons)-1):
        missense_prob = 0
        codon = codons[c]
        trinucl0 = f"{codons[c-1][2]}{codons[c][0:2]}"
        trinucl1 = f"{codons[c]}"
        trinucl2 = f"{codons[c][1:3]}{codons[c+1][0]}"

        # Print codon info
        if v:
            print(f"\n{colored('>>>')} Codon n: {c}", "\tRef AA:", aa, "\tCodon:", colored(codon))
            print("")
            print(f"\t\t{codons[c-1]}{colored(codons[c])}{codons[c+1]}")
            print(f"\t\t..{colored(trinucl0, 0, 0, 255)}....")
            print(f"\t\t...{colored(trinucl1, 0, 0, 255)}...")
            print(f"\t\t....{colored(trinucl2, 0, 0, 255)}..")
            print("")

        # Iterate through the possible contexts of a missense mut
        for i, trinucl in enumerate([trinucl0, trinucl1, trinucl2]):
            ref = trinucl[1]
            aa = gencode[codon]
            if v: 
                print(">> Context:", colored(trinucl, 0, 0, 255), "\n   Ref:", ref, "\n")
                print("Mut      Mut_prob                Alt_codon   Alt_AA")

            # Iterate through the possible alt 
            for alt in "ACGT":
                if alt != ref:         
                    alt_codon = [n for n in codon]
                    alt_codon[i] = alt
                    alt_codon = "".join(alt_codon)
                    alt_aa = gencode[alt_codon]  
                    # If there is a missense mut, get prob from context and sum it
                    if alt_aa != aa and alt_aa != "_":
                        mut = f"{trinucl}>{alt}"
                        if v: print(mut, "\t", mut_rate_dict[mut], "\t", alt_codon, "\t    ", alt_aa, )
                        missense_prob += mut_rate_dict[mut]
            if v: print("")

        if v: print(f">> Prob of missense mut: {missense_prob:.3}\n")
        missense_prob_vec.append(missense_prob)
    missense_prob_vec.append(0)                               # Assign 0 prob to the last residu

    # Convert into probabilities
    if get_probability:
        missense_prob_vec = np.array(missense_prob_vec) / sum(missense_prob_vec)
    
    return list(missense_prob_vec)


def get_miss_mut_prob_dict(mut_rate_dict, seq_df, v=False):
    """
    Given a dictionary of mut rate in 96 contexts (mut profile) and a 
    dataframe including Uniprot ID, HUGO symbol and DNA sequences, 
    get a dictionary with UniprotID-Fragment as keys and corresponding 
    vectors of missense mutation probabilities as values.
    """

    miss_prob_dict = {}
    # Process any Protein/fragment in the sequence df
    for _, row in seq_df.iterrows():
        if "F" in seq_df.columns:
            miss_prob_dict[f"{row.Uniprot_ID}-F{row.F}"] = get_miss_mut_prob(row.Seq_dna, mut_rate_dict, v=v)
        else:
            miss_prob_dict[f"{row.Uniprot_ID}"] = get_miss_mut_prob(row.Seq_dna, mut_rate_dict, v=v)

    return miss_prob_dict


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


    print("Input directory:", input_path)
    print("Output directory:", output_path)
    print("Path to DNA sequences df", path_seq_df, "\n")

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
        print(f" Processing {len(seq_df)} proteins/fragments in {filename}..")
        miss_prob_dict = get_miss_mut_prob_dict(mut_rate_dict=mut_rate_dict, seq_df=seq_df, v=verbose)

        # Save
        with open(f"{output_path}/{filename}.miss_mut_prob.json", "w") as fp:
            json.dump(miss_prob_dict, fp)


if __name__ == "__main__":
    main()
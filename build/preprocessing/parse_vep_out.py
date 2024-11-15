""" 
Run VEP and parse output for the 3D-clustering method.
"""


import pandas as pd
import numpy as np
import argparse


# Parser
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input_vcf", help = "Path of the vcf file used as input", type=str, required=True)
parser.add_argument("-o", "--output_path", help = "Path of the output file", type=str, required=True)
parser.add_argument("-c", "--canonical_only", help = "Include only canonical transcripts", type=bool, default=False)
args = parser.parse_args()
input_vcf = args.input_vcf
output_path = args.output_path
canonical_only = args.canonical_only


def read_vep_output(file_path):
    """
    Parse VEP output skipping heading comments.
    """
    
    with open(file_path, 'r') as file:
        # Skip lines starting with '##'
        lines = file.readlines()
        header_line = None
        for i, line in enumerate(lines):
            if not line.startswith("## "):
                header_line = i
                break
        if header_line is None:
            raise ValueError("No header line found.")

    vep_cols = ["#Uploaded_variation", "Feature", "Consequence", "SYMBOL", "CANONICAL", "Protein_position", "Amino_acids"]
    
    return pd.read_table(file_path, skiprows=header_line, usecols=vep_cols)


    cols_to_read = ["Variant_Classification",
                    "Tumor_Sample_Barcode",
                    "Feature", 
                    "Transcript_ID",
                    "Consequence", 
                    "SYMBOL", 
                    "Hugo_Symbol",
                    "CANONICAL", 
                    "HGVSp_Short",
                    "Amino_acids", 
                    "Protein_position"]


# Load vcf and rename columns for the output file
out_vep = read_vep_output(input_vcf)
out_vep = out_vep.rename(columns={"#Uploaded_variation" : "Tumor_Sample_Barcode",
                                  "Consequence" : "Variant_Classification",
                                  "SYMBOL" : "Hugo_Symbol",
                                  "Feature" : "Transcript_ID"})

# Select only mutations mapped to canonical transcript 
if canonical_only:
    out_vep = out_vep[out_vep["CANONICAL"] == "YES"]

# Parse mutation on the protein
out_vep["HGVSp_Short"] =  out_vep.apply(lambda x: 
                              "p." + x["Amino_acids"].split("/")[0] + str(x["Protein_position"]) + x["Amino_acids"].split("/")[1] 
                              if len(x["Amino_acids"].split("/")) > 1 
                              else np.nan, axis=1)
out_vep = out_vep[["Hugo_Symbol", "Tumor_Sample_Barcode", "Variant_Classification", "Transcript_ID", "HGVSp_Short"]]
out_vep.to_csv(output_path, sep="\t", index=False)

import re

import daiquiri
import numpy as np
import pandas as pd
import subprocess
import io
import gzip
import gc
import sys

from scripts import __logger_name__

logger = daiquiri.getLogger(__logger_name__ + ".run.utils")


## Parsers

def has_comments_as_header(filename):
    """
    Check if the file start with comments as headers.
    """
    
    if filename.endswith('.gz') or filename.endswith('.gzip'):
        with gzip.open(filename, 'rt') as file:
            for line in file:
                if line.startswith("## "):
                    return True
                else:
                    return False  
    else:
        with open(filename, 'r') as file:
            for line in file:
                if line.startswith("## "):
                    return True
                else:
                    return False
            
            
def read_csv_without_comments(path):
    """
    Read a csv file without any lines starting with '##' which 
    are commonly used as comments (e.g., as in the output of VEP).
    """
    
    # Run bash command
    command = f"cat {path} | grep -v '^## '"
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE)

    # Parse result
    captured_output = result.stdout.decode('utf-8')
    df = pd.read_csv(io.StringIO(captured_output), sep='\t', dtype={'Chromosome': str})
    
    return df


def get_seq_df_input_symbols(input_df, seq_df, mane=False):
    """
    Update gene names (HUGO Symbols) of O3D built sequence with names in input file.
    Do it only for entries in the sequence df with available transcript information
    and use transcript ID to get gene name.
    """

    # Split sequence df by entries with available transcript info (Reference_info 0 and 1) and not available ones (-1)
    seq_df_tr_missing = seq_df[seq_df["Reference_info"] == -1].reset_index(drop=True)
    seq_df_tr_available = seq_df[seq_df["Reference_info"] != -1].reset_index(drop=True)

    # Use names from input
    df_mapping = input_df[["Hugo_Symbol", "Feature"]].rename(columns={"Hugo_Symbol" : "Gene", "Feature" : "Ens_Transcr_ID"})
    seq_df_tr_available = seq_df_tr_available.drop(columns=["Gene"]).drop_duplicates().merge(df_mapping, how="left", on="Ens_Transcr_ID")

    # If the same gene is associated to multiple structures, keep the first one obtained from Uniprot (descending, Reference_info 1) or keep the MANE (ascending, Reference_info 0)
    # TO DO: Use the one reviewed (UniProtKB reviewed (Swiss-Prot)), if multiple Uniprot ones are present. The info must be added during the build step
    order_ascending = [True, mane]
    seq_df_tr_available = seq_df_tr_available.sort_values(by=["Gene", "Reference_info"], ascending=order_ascending).drop_duplicates(subset="Gene")

    # If the same genes is associated to multiple structures, keep the one not obtained by Backtranseq (Reference_info 1 or 0)
    seq_df = pd.concat([seq_df_tr_missing, seq_df_tr_available]).sort_values(by=["Gene", "Reference_info"], ascending=[True, False])

    return seq_df.drop_duplicates(subset="Gene").reset_index(drop=True)


def get_hgvsp_mut(df_row):
    """
    Parse mutation entries to get HGVSp_Short format.
    """
    
    amino_acids = df_row["Amino_acids"]
    
    if pd.isna(amino_acids):
        return np.nan
    
    amino_acids = amino_acids.split("/")
    if len(amino_acids) > 1:
        return f"p.{amino_acids[0]}{df_row['Protein_position']}{amino_acids[1]}"
    
    return np.nan
                
                
def filter_transcripts(df, seq_df):
    """
    Filter VEP output by Oncodrive3D transcripts. For genes with NA 
    transcripts in the sequence dataframe, keep canonical ones.
    """
    
    if "CANONICAL" in df.columns and "Feature" in df.columns:
        
        # Genes without available transcript info in O3D built datasets
        df_tr_missing = df[df["Hugo_Symbol"].isin(seq_df.loc[seq_df["Reference_info"] == -1, "Gene"])]
        df_tr_missing = df_tr_missing[df_tr_missing["CANONICAL"] == "YES"]
        
        # Genes with transcript info
        df_tr_available = df[df["Feature"].isin(seq_df.loc[seq_df["Reference_info"] != -1, "Ens_Transcr_ID"])]
        
        return pd.concat((df_tr_available, df_tr_missing))
    
    else:
        logger.critical("Failed to filter input by O3D transcripts. Please provide as input the output of VEP with canonical and transcripts information: Exiting..")
        sys.exit(1)
    

def reduce_sample_id(identifier):
    """
    Extract the middle section of the tumour sample identifier. 
    """
    
    parts = identifier.split('_')
    if len(parts) > 2:
        return parts[2]
    else:
        return identifier
    

def parse_vep_output(df,
                     seq_df=None, 
                     use_o3d_transcripts=False, 
                     use_input_symbols=False, 
                     mane=False):
    """
    Parse the dataframe in case it is the direct output of VEP without any 
    processing. Rename the columns to match the fields name of a MAF file, 
    and select the canonical transcripts if multiple ones are present.
    """

    # TO DO: check for mut ID instead of tumor sample (IF VEP OUTPUT AS INPUT)

    df.rename(columns={"SYMBOL": "Hugo_Symbol",
                       "Consequence": "Variant_Classification"}, inplace=True)
            
    # Adapt HUGO_Symbol in seq_df to input file
    if seq_df is not None and use_input_symbols:
        logger.debug("Adapting Oncodrive3D HUGO Symbols of built datasets to input file..")
        seq_df = get_seq_df_input_symbols(df, seq_df, mane)
    
    # Transcripts filtering
    if use_o3d_transcripts and seq_df is not None:
        logger.debug("Filtering input by Oncodrive3D built transcripts..")
        df = filter_transcripts(df, seq_df)
    elif "CANONICAL" in df.columns:
        df = df[df["CANONICAL"] == "YES"]
            
    # Get HGVSp
    if "HGVSp_Short" not in df.columns and "Amino_acids" in df.columns and "Protein_position" in df.columns:
        df["HGVSp_Short"] = df.apply(get_hgvsp_mut, axis=1)
    
    if "Tumor_Sample_Barcode" in df.columns:
        df["Tumor_Sample_Barcode"] = df['Tumor_Sample_Barcode'].apply(reduce_sample_id)
        
    return df, seq_df


def parse_mutations(maf):
    """
    Parse HGVSp_Short in maf.
    """
    
    # Ensure the required 'HGVSp_Short' column is present and not empty
    if 'HGVSp_Short' not in maf.columns or maf['HGVSp_Short'].isnull().all():
        logger.critical("Missing or empty 'HGVSp_Short' column in input MAF data.")
        sys.exit(1)

    # Parse the position, wild type, and mutation type from 'HGVSp_Short'
    maf.dropna(subset="HGVSp_Short", inplace=True)
    maf['Pos'] = maf['HGVSp_Short'].apply(lambda x: re.sub(r"\D", "", x)).astype(np.int32)
    maf['WT'] = maf['HGVSp_Short'].apply(lambda x: re.findall(r"\D", x)[2])
    maf['Mut'] = maf['HGVSp_Short'].apply(lambda x: re.findall(r"\D", x)[3])

    # Parse cols
    columns_to_keep = ['Hugo_Symbol', 'Pos', 'WT', 'Mut', 'Tumor_Sample_Barcode', 'Feature', 'Transcript_ID']
    columns_to_keep = [col for col in columns_to_keep if col in maf.columns]
    maf = maf[columns_to_keep]
    maf = maf.rename(columns={'Hugo_Symbol' : 'Gene', 'Feature': 'Transcript_ID'})

    return maf.sort_values(by=['Gene', 'Pos']).reset_index(drop=True)


def add_transcript_info(maf, seq_df):
    """
    Add transcript status information.
    """
    
    if 'Transcript_ID' not in maf.columns:
        maf['Transcript_ID'] = np.nan
    maf = maf.merge(seq_df[[col for col in ['Gene', 'Ens_Transcr_ID', 'Refseq_prot'] if col in seq_df.columns]].drop_duplicates(), 
                    on='Gene', how='left').rename(columns={"Ens_Transcr_ID" : "O3D_transcript_ID"})

    # Vectorized conditions for setting Transcript_status
    conditions = [
        maf['Transcript_ID'].isna(),
        maf['O3D_transcript_ID'].isna(),
        maf['Transcript_ID'] != maf['O3D_transcript_ID'],
        maf['Transcript_ID'] == maf['O3D_transcript_ID']
    ]
    choices = ['Input_missing', 'O3D_missing', 'Mismatch', 'Match']
    maf['Transcript_status'] = np.select(conditions, choices, default=np.nan)

    # Log transcript report
    transcript_report = maf['Transcript_status'].value_counts().reset_index(name='Count')
    transcript_report = ", ".join([f"{status} = {count}" for status, count in transcript_report.to_numpy()])
    logger.info(f"Transcript status of {len(maf)} mutations: {transcript_report}")

    return maf


def read_input(input_path):
    """
    Read input file optimizing memory usage.
    """

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
    
    header = pd.read_csv(input_path, delimiter='\t', nrows=0)
    cols_to_read = [col for col in cols_to_read if col in header.columns]
    dtype_mapping = {col : "object" for col in cols_to_read}

    # Filter the dtype mapping to include only the columns you are going to read
    dtype = {key: dtype_mapping[key] for key in cols_to_read if key in dtype_mapping}
    
    return pd.read_csv(input_path, delimiter='\t', usecols=cols_to_read, dtype=dtype)


def parse_maf_input(maf_input_path, seq_df=None, use_o3d_transcripts=False, use_input_symbols=False, mane=False):
    """
    Parsing and process MAF input data.
    """

    # Load, parse from VEP and update seq_df if needed
    logger.info(f"Reading input mutations file..")
    maf = read_input(maf_input_path)
    logger.debug(f"Processing [{len(maf)}] total mutations..")
    maf, seq_df = parse_vep_output(maf, seq_df, use_o3d_transcripts, use_input_symbols, mane)
    
    # ### DEBUG COMPARE INPUT
    # logger.critical("Writing maf_00.tsv for debugging")
    # maf.to_csv("/workspace/projects/clustering_3d/o3d_analysys/datasets/output/normal/o3d_output/2024/maf_00.tsv", sep="\t")
    
    # ###

    # Extract and parse missense mutations
    maf = maf[maf['Variant_Classification'].str.contains('Missense_Mutation|missense_variant')]
    if "Protein_position" in maf.columns:
        maf = maf[~maf['Protein_position'].astype(str).str.contains('-')] # Filter DBS
    logger.debug(f"Processing [{len(maf)}] missense mutations..")
    maf = parse_mutations(maf) 
    
    # Add transcript status from seq_df
    if seq_df is not None:
        maf = add_transcript_info(maf, seq_df)

    # ### DEBUG COMPARE INPUT
    # logger.critical("Writing maf_01.tsv for debugging")
    # maf.to_csv("/workspace/projects/clustering_3d/o3d_analysys/datasets/output/normal/o3d_output/2024/maf_01.tsv", sep="\t")
    # ###
    
    return maf.reset_index(drop=True), seq_df


## Other utils

def get_gene_entry(data, genes, entry):

    return [data.loc[data["Gene"] == gene, entry].values[0] for gene in genes]


def weighted_avg_plddt_vol(target_pos, mut_plddt_df, cmap):
    """
    Get the weighted average pLDDT score across residues in 
    the volume, based on number of mutations hitting each residue.  
    """
    
    return mut_plddt_df[[pos in np.where(cmap[target_pos-1])[0]+1 for pos in mut_plddt_df.Pos]].Confidence.mean()


def weighted_avg_pae_vol(target_pos, mut_plddt_df, cmap, pae):
    """
    Get the weighted average PAE across residues in the volume, 
    based on number of mutations hitting each residue.  
    """
    
    contacts = mut_plddt_df[[pos in np.where(cmap[target_pos-1])[0]+1 for pos in mut_plddt_df.Pos]].Pos.values
    pae_vol = pae[np.repeat(target_pos-1, len(contacts)), contacts-1].mean()
    
    return pae_vol


def get_samples_info(mut_gene_df, cmap):
    """
    Get total samples and ratio of unique samples having
    mutations in the volume of each mutated residues.
    """
    
    # Get total samples and # mutated samples of each mutated res
    tot_samples = len(mut_gene_df["Tumor_Sample_Barcode"].unique())
    pos_barcodes = mut_gene_df.groupby("Pos").apply(lambda x: x["Tumor_Sample_Barcode"].unique())
    pos_barcodes = pos_barcodes.reset_index().rename(columns={0 : "Barcode"})

    # Get the ratio of unique samples with mut in the vol of each mutated res
    uniq_pos_barcodes = [len(pos_barcodes[[pos in np.where(cmap[i-1])[0]+1 for 
                                        pos in pos_barcodes.Pos]].Barcode.explode().unique()) for i in pos_barcodes.Pos]
    pos_barcodes["Tot_samples"] = tot_samples        
    pos_barcodes["Samples_in_vol"] = uniq_pos_barcodes
    #pos_barcodes["Ratio_samples_in_vol"] = np.array(uniq_pos_barcodes) / tot_samples     
    
    return pos_barcodes


# def get_unique_pos_in_contact(lst_pos, cmap):                                                     # ----> TO TEST <----
#     """
#     Identify unique positions in contact with given ones.
#     """
#     # Adjust indices for 0-based indexing
#     indices = np.array(lst_pos) - 1
    
#     # Use advanced indexing to extract rows from cmap based on indices
#     contact_positions = cmap[indices, :]
    
#     # Find all positions where contact exists, adjust for 1-based indexing
#     return np.unique(np.where(contact_positions)[1] + 1)


# def add_samples_info(mut_gene_df, result_pos_df, samples_info, cmap, pae=None):                   # ----> TO TEST <----
#     """
#     Add information about the ratio of unique samples in the volume of 
#     each mutated residue and in each detected community (meta-cluster).
#     """
    
#     # Merge samples information by position
#     result_pos_df = result_pos_df.merge(samples_info.drop(columns=["Barcode"]), on="Pos", how="outer")
    
#     if result_pos_df["Cluster"].notna().any():
#         # Compute unique positions in contact for each cluster
#         community_pos = result_pos_df.groupby("Cluster")["Pos"].agg(list)
#         community_contact_pos = community_pos.apply(lambda positions: get_unique_pos_in_contact(positions, cmap))
        
#         # Compute metrics for each community
#         def community_metrics(positions):
#             in_contact = mut_gene_df['Pos'].isin(positions)
#             samples = mut_gene_df.loc[in_contact, 'Tumor_Sample_Barcode'].unique()
#             mutations = in_contact.sum()
#             confidence = mut_gene_df.loc[in_contact, 'Confidence'].mean() if pae else np.nan
            
#             return pd.Series({
#                 "Mut_in_cl_vol": mutations,
#                 "Samples_in_cl_vol": len(samples),
#                 "Res_in_cl": len(positions),
#                 "pLDDT_cl_vol": confidence
#             })
        
#         metrics = community_contact_pos.apply(community_metrics)
#         result_pos_df = result_pos_df.join(metrics, on="Cluster", how="left")
#     else:
#         result_pos_df[["Samples_in_cl_vol", "Mut_in_cl_vol", "Res_in_cl", "pLDDT_cl_vol"]] = np.nan
    
#     return result_pos_df


def get_unique_pos_in_contact(lst_pos, cmap):
    """
    Given a list of position and a contact map, return a numpy 
    array of unique positions in contact with the given ones.
    """
    
    return np.unique(np.concatenate([np.where(cmap[pos-1])[0]+1 for pos in lst_pos]))


def add_info(mut_gene_df, result_pos_df, cmap, pae=None):
    """
    Add information about the ratio of unique samples in the volume of 
    each mutated residues and in each detected community (meta-cluster) 
    to the residues-level output of the tool.
    """

    # Add sample info
    if "Tumor_Sample_Barcode" in mut_gene_df.columns:
        samples_info = get_samples_info(mut_gene_df, cmap)
        result_pos_df = result_pos_df.merge(samples_info.drop(columns=["Barcode"]), on="Pos", how="outer")
    else:
        result_pos_df["Tot_samples"] = np.nan
        result_pos_df["Samples_in_vol"] = np.nan

    # Get per-community info
    if result_pos_df["Cluster"].isna().all():
        result_pos_df["Samples_in_cl_vol"] = np.nan
        result_pos_df["Mut_in_cl_vol"] = np.nan
        result_pos_df["Res_in_cl"] = np.nan
        result_pos_df["pLDDT_cl_vol"] = np.nan
    else:       
        community_pos = result_pos_df.groupby("Cluster").apply(lambda x: x.Pos.values)
        community_mut = community_pos.apply(lambda x: sum([pos in get_unique_pos_in_contact(x, cmap) for 
                                                           pos in mut_gene_df.Pos]))
        community_plddt = community_pos.apply(lambda x: mut_gene_df.Confidence[[pos in get_unique_pos_in_contact(x, cmap) 
                                                                             for pos in mut_gene_df.Pos]].mean())
        community_pos_count = community_pos.apply(lambda x: len(x))
        if "Tumor_Sample_Barcode" in mut_gene_df.columns:
            community_samples = community_pos.apply(lambda x: 
                                            len(mut_gene_df[[pos in get_unique_pos_in_contact(x, cmap) for 
                                                            pos in mut_gene_df.Pos]].Tumor_Sample_Barcode.unique()))
        else:
            community_samples = np.nan
        community_samples = pd.DataFrame({"Samples_in_cl_vol" : community_samples, 
                                          "Mut_in_cl_vol" : community_mut,
                                          "Res_in_cl" : community_pos_count,
                                          "pLDDT_cl_vol" : community_plddt})
        
        # Add to residues-level result
        result_pos_df = result_pos_df.merge(community_samples, on="Cluster", how="outer")
        
    # AF PAE
    if pae is not None:
        result_pos_df["PAE_vol"] = result_pos_df.apply(lambda x: weighted_avg_pae_vol(x["Pos"], mut_gene_df, cmap, pae), axis=1)  
    else:
        result_pos_df["PAE_vol"] = np.nan
        
    # AF confidence
    result_pos_df["pLDDT_res"] = result_pos_df.apply(lambda x: mut_gene_df.Confidence[mut_gene_df["Pos"] == x.Pos].values[0], axis=1)
    result_pos_df["pLDDT_vol"] = result_pos_df.apply(lambda x: weighted_avg_plddt_vol(x["Pos"], mut_gene_df, cmap), axis=1)
    result_pos_df["pLDDT_cl_vol"] = result_pos_df.pop("pLDDT_cl_vol")

    # Sort positions
    result_pos_df = result_pos_df.sort_values("Rank").reset_index(drop=True)
    
    return result_pos_df


def add_nan_clust_cols(result_gene):
    """
    Add columns showing clustering results with only NA for 
    genes that are not tested (not enough mutations, etc).
    """

    result_gene = result_gene.copy()
    
    columns = ["pval", 
               "qval", 
               "C_gene", 
               "C_pos", 
               'C_label', 
               'Ratio_obs_sim_top_vol', 
               "Clust_res", 
               'Clust_mut', 
               'Mut_in_top_vol', 
               "Mut_in_top_cl_vol",
               'Tot_samples', 
               'Samples_in_top_vol', 
               'Samples_in_top_cl_vol', 
               "PAE_top_vol", 
               "pLDDT_top_vol", 
               "pLDDT_top_cl_vol", 
               'F']
    
    for col in columns:
        result_gene[col] = np.nan

    return result_gene


def sort_cols(result_gene):
    """
    Simply change the order of columns of the genes-level result.
    """

    cols = ['Gene', 
            'HGNC_ID',
            'Ens_Gene_ID', 
            'Refseq_prot',
            'Uniprot_ID', 
            'pval', 
            'qval', 
            'C_gene', 
            'C_pos', 
            'C_label', 
            'Ratio_obs_sim_top_vol', 
            "Clust_res",
            'Mut_in_gene', 
            'Clust_mut', 
            'Mut_in_top_vol', 
            "Mut_in_top_cl_vol",
            'Tot_samples', 
            'Samples_in_top_vol', 
            'Samples_in_top_cl_vol', 
            "PAE_top_vol", 
            "pLDDT_top_vol", 
            "pLDDT_top_cl_vol",
            'F', 
            'Ratio_not_in_structure',
            'Ratio_WT_mismatch',
            'Mut_zero_mut_prob',
            'Pos_zero_mut_prob',
            'Ratio_mut_zero_prob',
            'Cancer', 
            'Cohort',
            'Transcript_ID',
            'O3D_transcript_ID',
            'Transcript_status',
            'Status']

    return result_gene[[col for col in cols if col in result_gene.columns]]


def empty_result_pos():
    """
    Get an empty position-level result of the clustering method.
    """
    
    cols = ['Gene', 
            'Uniprot_ID', 
            'F', 
            'Pos', 
            'Mut_in_gene', 
            'Mut_in_res',
            'Mut_in_vol', 
            'Obs_anomaly', 
            'Ratio_obs_sim', 
            'pval', 
            'C', 
            'C_ext',
            'Cluster', 
            'Rank', 
            'Tot_samples', 
            'Samples_in_vol', 
            'Samples_in_cl_vol',
            'Mut_in_cl_vol', 
            'Res_in_cl', 
            'PAE_vol', 
            'pLDDT_res', 
            'pLDDT_vol',
            'pLDDT_cl_vol', 
            'Cancer', 
            'Cohort']
    
    return pd.DataFrame(columns=cols)
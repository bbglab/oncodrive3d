import re

import daiquiri
import numpy as np
import pandas as pd
import subprocess
import io
import gzip

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
    seq_df = seq_df.copy()
    seq_df_tr_missing = seq_df[seq_df["Reference_info"] == -1].reset_index(drop=True)
    seq_df_tr_available = seq_df[seq_df["Reference_info"] != -1].reset_index(drop=True)

    # Use names from input
    seq_df_tr_available = seq_df_tr_available.drop(columns=["Gene"]).drop_duplicates()
    df_mapping = input_df[["Hugo_Symbol", "Feature"]].rename(columns={"Hugo_Symbol" : "Gene", "Feature" : "Ens_Transcr_ID"})
    seq_df_tr_available = seq_df_tr_available.merge(df_mapping, how="left", on="Ens_Transcr_ID")

    # If the same gene is associated to multiple structures, keep the first one obtained from Uniprot (descending, Reference_info 1) or keep the MANE (ascending, Reference_info 0)
    # TO DO: Use the one reviewed (UniProtKB reviewed (Swiss-Prot)), if multiple Uniprot ones are present. The info must be added during the build step
    if mane:
        ascending = True
    else:
        ascending = False
    seq_df_tr_available = seq_df_tr_available.sort_values(["Gene", "Reference_info"], ascending=ascending).drop_duplicates(subset='Gene').reset_index(drop=True)
    
    seq_df = pd.concat((seq_df_tr_missing, seq_df_tr_available)).sort_values(["Gene"]).reset_index(drop=True)

    # If the same genes is associated to multiple structures, keep the one not obtained by Backtranseq (Reference_info -1)
    seq_df = seq_df.sort_values(["Gene", "Reference_info"], ascending=ascending).drop_duplicates(subset='Gene').reset_index(drop=True)

    return seq_df


def get_hgvsp_mut(df_row):
    """
    Parse mutation entries to get HGVSp_Short format.
    """
    
    if pd.isna(df_row["Amino_acids"]):
        return np.nan
    
    elif len(df_row["Amino_acids"].split("/")) > 1:
        return "p." + df_row["Amino_acids"].split("/")[0] + str(df_row["Protein_position"]) + df_row["Amino_acids"].split("/")[1] 
    
    else:
        return np.nan


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
    
    rename_dict = {"SYMBOL": "Hugo_Symbol", 
                   "Consequence": "Variant_Classification",
                   "#Uploaded_variation" : "Tumor_Sample_Barcode",
                   "#UploadedVariation" : "Tumor_Sample_Barcode"}

    for key, value in rename_dict.items():
        if key in df.columns and value not in df.columns:
            df.rename(columns={key: value}, inplace=True)
            
    # Adapt Hugo_Symbol in seq_df to input file
    if seq_df is not None and use_input_symbols:
        logger.debug("Adapting Oncodrive3D HUGO Symbols of built datasets to input file..")
        seq_df = get_seq_df_input_symbols(df, seq_df, mane)
    
    ## Transcripts filtering
    
    # Include O3D transcripts
    if use_o3d_transcripts:
        logger.info("Filtering input by Oncodrive3D built transcripts..")
        if seq_df is not None:
            if "CANONICAL" in df.columns and "Feature" in df.columns and "Transcript" in df["Feature_type"].unique():
                # For the genes without available transcript info in O3D built datasets, select canonical
                df_tr_missing = df[df["Hugo_Symbol"].isin(seq_df.loc[seq_df["Reference_info"] == -1, "Gene"])]
                df_tr_missing = df_tr_missing[df_tr_missing["CANONICAL"] == "YES"]
                # For those with transcript info, select transcript in O3D build datasets
                df_tr_available = df[df["Feature"].isin(seq_df.loc[seq_df["Reference_info"] != -1, "Ens_Transcr_ID"])]
                df = pd.concat((df_tr_available, df_tr_missing))
            else:
                raise RuntimeError(f"Failed to filter input by O3D transcripts. Please provide as input the output of VEP with canonical and transcripts information: Exiting..")
        else:
            raise RuntimeError(f"Failed to filter input by O3D transcripts. Dataframe of sequences not provided: Exiting..")
    else:
        # Include canonical transcripts
        if "CANONICAL" in df.columns:
            if (df["CANONICAL"] != "YES").any():
                logger.warning("The 'CANONICAL' field is present in the MAF file: Selecting only canonical transcripts..")
                df = df[df["CANONICAL"] == "YES"]
            
    # Get HGVSp
    if "HGVSp_Short" not in df.columns:
        if "Amino_acids" in df.columns and "Protein_position" in df.columns:
            logger.debug("Input detected as direct VEP output: Parsing translational effect of variant..")            
            df["HGVSp_Short"] = df.apply(lambda x: get_hgvsp_mut(x), axis=1)
        else:
            logger.critical("Translational effect of variant allele not found: Input file must include either the field 'HGVSp_Short' or 'Amino_acids' and 'Protein_positions'.")
        
    return df, seq_df


def parse_maf_input(maf_input_path, 
                    seq_df=None, 
                    use_o3d_transcripts=False, 
                    use_input_symbols=False, 
                    mane=False):
    """
    Parse the MAF file which is used as input for Oncodrive3D.
    """

    # Load, parse from VEP and update seq_df if needed
    logger.debug(f"Reading input MAF...")
    if has_comments_as_header(maf_input_path):
        logger.debug("Detected '##' comments in the file header: Reading the file without comments...")
        maf = read_csv_without_comments(maf_input_path)
    else:
        maf = pd.read_csv(maf_input_path, sep="\t", dtype={'Chromosome': str})
    logger.debug(f"Processing [{len(maf)}] total mutations..")
    maf, seq_df = parse_vep_output(maf, 
                                   seq_df, 
                                   use_o3d_transcripts, 
                                   use_input_symbols, 
                                   mane)

    # Select only missense mutation and extract mutations
    maf = maf[maf['Variant_Classification'].str.contains('Missense_Mutation')
              | maf['Variant_Classification'].str.contains('missense_variant')]
    
    # TODO: maybe change the parsing process and make it simpler
    # TODO: DBS filtering only occurs if the "Protein_position" col is in the input
    #       maybe change such that it filter them out anyway
    
    # Filtering DBS                                
    if "Protein_position" in maf.columns:
        dbs_ix = maf.Protein_position.apply(lambda x: len(str(x).split("-"))) > 1
        if sum(dbs_ix) > 0:
            logger.warning(f"Parsed {sum(dbs_ix)} double base substitutions: Filtering the mutations...")
            maf = maf[~dbs_ix]
    
    # Parse mut
    logger.debug(f"Processing [{len(maf)}] missense mutations...")
    maf = maf.dropna(subset="HGVSp_Short")
    maf["Pos"] = maf.loc[:, "HGVSp_Short"].apply(lambda x: int(re.sub("\\D", "", (x[2:]))))
    maf["WT"] = maf["HGVSp_Short"].apply(lambda x: re.findall("\\D", x[2:])[0])
    maf["Mut"] = maf["HGVSp_Short"].apply(lambda x: re.findall("\\D", x[2:])[1])
    maf = maf[["Hugo_Symbol", "Pos", "WT", "Mut"] + [col for col in ["Tumor_Sample_Barcode", "Feature", "Transcript_ID"] if col in maf.columns]]
    maf = maf.sort_values("Pos").rename(columns={"Hugo_Symbol" : "Gene"})
    if "Feature" in maf.columns:
        maf = maf.rename(columns={"Feature" : "Transcript_ID"})
    
    # Add transcript info
    if seq_df is not None:
        if "Transcript_ID" not in maf.columns:
            logger.warning("Input transcript ID not found")
            maf["Transcript_ID"] = np.nan
        maf = maf.merge(seq_df[[col for col in ["Gene", "Ens_Transcr_ID", "Refseq_prot"] if col in seq_df.columns]], how="left", on="Gene")
        maf = maf.rename(columns={"Ens_Transcr_ID" : "O3D_transcript_ID"})
        maf["Transcript_status"] = maf.apply(lambda x: "Input_missing" if pd.isna(x["Transcript_ID"]) 
                                             else "O3D_missing" if pd.isna(x["O3D_transcript_ID"]) 
                                             else "Mismatch" if x["Transcript_ID"] != x["O3D_transcript_ID"]
                                             else "Match" if x["Transcript_ID"] == x["O3D_transcript_ID"]
                                             else np.nan, axis=1)
        
        transcript_report = maf.Transcript_status.value_counts().reset_index()
        transcript_report.columns = "Status", "Count"
        transcript_report = ", ".join([f"{v}: {c}" for (v, c) in zip(transcript_report.Status, 
                                                                     transcript_report.Count)])
        logger.info(f"Tanscript status = {transcript_report}")
    
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


def get_unique_pos_in_contact(lst_pos, cmap):
    """
    Given a list of position and a contact map, return a numpy 
    array of unique positions in contact with the given ones.
    """
    
    return np.unique(np.concatenate([np.where(cmap[pos-1])[0]+1 for pos in lst_pos]))


def add_samples_info(mut_gene_df, result_pos_df, samples_info, cmap, pae=None):
    """
    Add information about the ratio of unique samples in the volume of 
    each mutated residues and in each detected community (meta-cluster) 
    to the residues-level output of the tool.
    """

    # Add samples in vol
    result_pos_df = result_pos_df.merge(samples_info.drop(columns=["Barcode"]), on="Pos", how="outer")
    
    # Get per-community ratio of mutated samples
    if result_pos_df["Cluster"].isna().all():
        result_pos_df["Samples_in_cl_vol"] = np.nan
        result_pos_df["Mut_in_cl_vol"] = np.nan
        result_pos_df["Res_in_cl"] = np.nan
        result_pos_df["pLDDT_cl_vol"] = np.nan
        #result_pos_df["Ratio_samples_in_comm"] = np.nan
    else:       
        community_pos = result_pos_df.groupby("Cluster").apply(lambda x: x.Pos.values)
        community_mut = community_pos.apply(lambda x: sum([pos in get_unique_pos_in_contact(x, cmap) for 
                                                           pos in mut_gene_df.Pos]))
        community_samples = community_pos.apply(lambda x: 
                                        len(mut_gene_df[[pos in get_unique_pos_in_contact(x, cmap) for 
                                                         pos in mut_gene_df.Pos]].Tumor_Sample_Barcode.unique()))
        community_plddt = community_pos.apply(lambda x: mut_gene_df.Confidence[[pos in get_unique_pos_in_contact(x, cmap) 
                                                                             for pos in mut_gene_df.Pos]].mean())

        #community_samples_ratio = community_samples / samples_info["Tot_samples"].unique()[0]
        #community_mut_ratio = community_mut / len(mut_gene_df)
        community_pos_count = community_pos.apply(lambda x: len(x))
        community_samples = pd.DataFrame({"Samples_in_cl_vol" : community_samples, 
                                          #"Ratio_samples_in_comm" : community_samples_ratio,
                                          "Mut_in_cl_vol" : community_mut,
                                          #"Ratio_mut_in_comm" : community_mut_ratio,
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
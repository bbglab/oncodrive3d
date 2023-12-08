import pandas as pd
from itertools import product

import daiquiri
import numpy as np
import json

from scripts.run.mutability import init_mutabilities_module
from scripts.run.miss_mut_prob import get_miss_mut_prob_dict, mut_rate_vec_to_dict, get_unif_gene_miss_prob
from scripts import __logger_name__

logger = daiquiri.getLogger(__logger_name__ + ".annotations.plot")



# ==============================
# Missense mutations probability
# ==============================

def get_miss_mut_prob(seq_df,
                      mutability_config_path, 
                      mut_profile_path):
    
    # Missense mut prob with mutabilities
    if mutability_config_path is not None:
        print(f"Computing missense mut probabilities using mutabilities...")
        mutab_config = json.load(open(mutability_config_path))
        init_mutabilities_module(mutab_config)
        seq_df = seq_df[seq_df["Reference_info"] == 1]   
        seq_df['Exons_coord'] = seq_df['Exons_coord'].apply(eval)  
        genes_to_process = [gene for gene in genes_to_process if gene in seq_df["Gene"].unique()]
        miss_prob_dict = get_miss_mut_prob_dict(mut_rate_dict=None, seq_df=seq_df,
                                                mutability=True, mutability_config=mutab_config)
        
    # Without mutabilities
    elif mut_profile_path is not None:
        # Compute dict from mut profile of the cohort and dna sequences
        mut_profile = json.load(open(mut_profile_path))
        print(f"Computing missense mut probabilities...")
        if not isinstance(mut_profile, dict):
            mut_profile = mut_rate_vec_to_dict(mut_profile)
        miss_prob_dict = get_miss_mut_prob_dict(mut_rate_dict=mut_profile, seq_df=seq_df)
    
    # Uniform distribution
    else:
        miss_prob_dict = None
        
        
# ============
# PLOT WRAPPER
# ============

def generate_plot(datasets):
    
    seq_df_path = f"{datasets}/seq_for_mut_prob.csv"
    seq_df = pd.read_csv(seq_df_path)
    
    
    ## Select the genes from gene_result
    # IF: processed only
    # IF: hits only
    # IF: top
    
    ## Subset seq_df
    
    ## Get miss mut prob 
    
    ## Generate directory output
    
    ## Summary plot
    
    ## Generate directory output/gene_plots
    
    ## Diagnostic plot
    
    pass
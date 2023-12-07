import pandas as pd



# ## Get missense mutations probability

# data = parse_maf_input(input_maf_path)
# if len(data) > 0:
#     seq_df = pd.read_csv(seq_df_path)
#     plddt_df = pd.read_csv(plddt_path, dtype={"Pos" : int, 
#                                             "Res" : str, 
#                                             "Confidence" : float, 
#                                             "Uniprot_ID" : str, 
#                                             "AF_F" : str})


#     ## Run

#     result_np_gene_lst = []

#     # Get genes with enough mut
#     genes = data.groupby("Gene").apply(lambda x: len(x))
#     genes_mut = genes[genes >= 2]
#     genes_no_mut = genes[genes < 2].index
#     if len(genes_no_mut) > 0:
#         print(f"Detected [{len(genes_no_mut)}] genes without enough mutations: Skipping...")
#         result_gene = pd.DataFrame({"Gene" : genes_no_mut,
#                                     "Uniprot_ID" : np.nan,
#                                     "F" : np.nan,
#                                     "Mut_in_gene" : 1,
#                                     "Max_mut_pos" : np.nan,
#                                     "Structure_max_pos" : np.nan,
#                                     "Status" : "No_mut"})
#         result_np_gene_lst.append(result_gene)   

#     # Get genes with corresponding Uniprot-ID mapping
#     gene_to_uniprot_dict = {gene : uni_id for gene, uni_id in seq_df[["Gene", "Uniprot_ID"]].drop_duplicates().values}
#     genes_to_process = [gene for gene in genes_mut.index if gene in gene_to_uniprot_dict.keys()]
#     seq_df = seq_df[[gene in genes_to_process for gene in seq_df["Gene"]]].reset_index(drop=True)
#     genes_no_mapping = genes[[gene in genes_mut.index and gene not in gene_to_uniprot_dict.keys() for gene in genes.index]]
#     if len(genes_no_mapping) > 0:
#         print(f"Detected [{len(genes_no_mapping)}] genes without IDs mapping: Skipping...")
#         result_gene = pd.DataFrame({"Gene" : genes_no_mapping.index,
#                                     "Uniprot_ID" : np.nan,
#                                     "F" : np.nan,
#                                     "Mut_in_gene" : genes_no_mapping.values,
#                                     "Max_mut_pos" : np.nan,
#                                     "Structure_max_pos" : np.nan,
#                                     "Status" : "No_ID_mapping"})
#         result_np_gene_lst.append(result_gene)
    
#     # Filter on fragmented (AF-F) genes
#     if no_fragments:
#         # Return the fragmented genes as non processed output
#         genes_frag = seq_df[seq_df.F.str.extract(r'(\d+)', expand=False).astype(int) > 1]
#         genes_frag = genes_frag.Gene.reset_index(drop=True).values
#         genes_frag_mut = genes_mut[[gene in genes_frag for gene in genes_mut.index]]
#         genes_frag = genes_frag_mut.index.values
#         genes_frag_id = [gene_to_uniprot_dict[gene] for gene in genes_frag]     
#         if len(genes_frag) > 0:   
#             print(f"Detected [{len(genes_frag)}] fragmented genes with disabled fragments processing: Skipping...")
#             result_gene = pd.DataFrame({"Gene" : genes_frag,
#                                         "Uniprot_ID" : genes_frag_id,
#                                         "F" : np.nan,
#                                         "Mut_in_gene" : genes_frag_mut.values,
#                                         "Max_mut_pos" : np.nan,
#                                         "Structure_max_pos" : np.nan,
#                                         "Status" : "Fragmented"})
#             result_np_gene_lst.append(result_gene)
#             # Filter out from genes to process and seq df
#             genes_to_process = [gene for gene in genes_to_process if gene not in genes_frag]
#             seq_df = seq_df[[gene in genes_to_process for gene in seq_df["Gene"]]].reset_index(drop=True)
            

#     # Missense mut prob
#     # using mutabilities if provided
#     if  mutability_config_path is not None:
#         print(f"Computing missense mut probabilities using mutabilities...")
#         mutab_config = json.load(open(mutability_config_path))
#         init_mutabilities_module(mutab_config)
#         seq_df = seq_df[seq_df["Reference_info"] == 1]   
#         seq_df['Exons_coord'] = seq_df['Exons_coord'].apply(eval)  
#         genes_to_process = [gene for gene in genes_to_process if gene in seq_df["Gene"].unique()]
#         genes_not_mutability = [gene for gene in genes_to_process if gene not in seq_df["Gene"].unique()]
#         miss_prob_dict = get_miss_mut_prob_dict(mut_rate_dict=None, seq_df=seq_df,
#                                                 mutability=True, mutability_config=mutab_config)
        
#         # TODO: return Uniprot_ID, F, etc
#         if len(genes_not_mutability) > 0:   
#             print(f"Detected [{len(genes_not_mutability)}] genes without mutability information: Skipping...")
#             result_gene = pd.DataFrame({"Gene" : genes_not_mutability,
#                                         "Uniprot_ID" : np.nan,
#                                         "F" : np.nan,
#                                         "Mut_in_gene" : np.nan,
#                                         "Max_mut_pos" : np.nan,
#                                         "Structure_max_pos" : np.nan,
#                                         "Status" : "No_mutability"})
#             result_np_gene_lst.append(result_gene)
#     # using mutational profiles
#     elif mut_profile_path is not None:
#         # Compute dict from mut profile of the cohort and dna sequences
#         mut_profile = json.load(open(mut_profile_path))
#         print(f"Computing missense mut probabilities...")
#         if not isinstance(mut_profile, dict):
#             mut_profile = mut_rate_vec_to_dict(mut_profile)
#         miss_prob_dict = get_miss_mut_prob_dict(mut_rate_dict=mut_profile, seq_df=seq_df)
#     else:
#         print(f"Mutation profile not provided: Uniform distribution will be used for scoring and simulations.")
#         miss_prob_dict = None

def plot():
    
    pass
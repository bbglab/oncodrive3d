"""
Use ChimeraX to generate 3D structures colored by metrics

singularity exec ~/Repositories/containers/chimerax/o3d_chimerax.sif python3 /o3d_chimerax_plot.py \
    -o /home/spellegrini/PhD/Projects/Oncodrive3D/normal_tissue/bladder \
        -g /home/spellegrini/PhD/Projects/Oncodrive3D/normal_tissue/bladder/BLCA_NORMAL.3d_clustering_genes.csv \
            -p /home/spellegrini/PhD/Projects/Oncodrive3D/normal_tissue/bladder/BLCA_NORMAL.3d_clustering_pos.csv \
                -d /home/spellegrini/PhD/Projects/Oncodrive3D/normal_tissue/bladder \
                    -s /home/spellegrini/PhD/Projects/Oncodrive3D/normal_tissue/bladder/BLCA_NORMAL.seq_df.processed.tsv \
                        -c BLCA_NORMAL \
                            --fragmented_proteins \
                                --transparent_bg
"""

import pandas as pd   
import numpy as np
import os 
import subprocess
import math
import click
import logging
import daiquiri

from scripts import __logger_name__
logger = daiquiri.getLogger(__logger_name__ + ".plotting.chimerax_plot")


def create_attribute_file(path_to_file, 
                        df, 
                        attribute_col,
                        pos_col="Pos",
                        attribute_name="local_attribute",
                        gene="Gene_name", 
                        protein="Protein_name"):

    logger.info(f"Saving {path_to_file}")
    with open(path_to_file, 'w') as f:
        f.write("#\n")
        f.write(f"#  Mutations in volume for {protein} ({gene})\n")
        f.write("#\n")
        f.write("#  Use this file to assign the attribute in Chimera with the \n")
        f.write("#  Define Attribute tool or the command defattr.\n")
        f.write("#\n")
        f.write(f"attribute: {attribute_name}\n")
        f.write("recipient: residues\n")

        for _, row in df.iterrows():
            f.write(f"\t:{int(row[pos_col])}\t{row[attribute_col]}\n")
            
            
def round_to_first_nonzero(num):
    if num == 0:
        return 0 
    
    scale = -int(math.floor(math.log10(abs(num))))
    shifted_num = num * (10 ** scale)
    rounded_num = round(shifted_num)
    result = rounded_num / (10 ** scale)
    
    return result
            

def is_float_an_integer(value):
    if isinstance(value, float):
        return value.is_integer()
    return False


def get_intervals(attribute_vector, attribute):
    
    max_value = attribute_vector.max()
    min_value = attribute_vector.min()
    min_value = 0 if attribute == "score" else 1
    intervals = np.linspace(min_value, max_value, 5)
    intervals = np.round(intervals, 2) if attribute == "score" else np.round(intervals)
    if attribute == "logscore":
        pos_values = np.linspace(0, max_value, 3)
        intervals = np.round([-max_value, -pos_values[1], 0, pos_values[1], max_value], 2)
    intervals = [round(n) if is_float_an_integer(n) else n for n in intervals]
        
    return intervals


def get_key_logscore(intervals):
    return f"{intervals[0]},#0571B0:{intervals[1]},#92C5DE:{intervals[2]},white:{intervals[3]},#F4A582:{intervals[4]},#CA0020"


def get_chimerax_command(chimerax_bin, 
                        pdb_path, 
                        chimera_output_path, 
                        attr_file_path, 
                        attribute, 
                        intervals, 
                        gene, 
                        uni_id,
                        labels,
                        i,
                        f,
                        cohort="",
                        clusters=None,
                        pixelsize=0.1,
                        palette="YlOrRd-5",
                        transparent_bg=False):
    
    chimerax_command = (
        f"{chimerax_bin} --nogui --offscreen --silent --cmd "
        f"\"open {pdb_path}; "
        "set bgColor white; "
        "color lightgray; "
        f"open {attr_file_path}; "
        f"color byattribute {attribute} palette {palette}; "
        "hide atoms;"
        "show cartoons;"
        f"key {palette} :{intervals[0]} :{intervals[1]}  :{intervals[2]} :{intervals[3]} :{intervals[4]} pos 0.35,0.03 fontSize 4 size 0.3,0.02;"
        f"2dlabels create label text '{labels[attribute]}' size 6 color darkred xpos 0.35 ypos 0.065;"
        f"2dlabels create title text '{gene} - {uni_id}-F{f} ' size 6 color darkred xpos 0.35 ypos 0.93;"
        "lighting soft;"
        "graphics silhouettes true width 1.3;"
        "zoom;"
    )
    
    if attribute == "logscore":
        chimerax_command = (
            f"{chimerax_bin} --nogui --offscreen --silent --cmd "
            f"\"open {pdb_path}; "
            "set bgColor white; "
            "color lightgray; "
            f"open {attr_file_path}; "
            f"color byattribute {attribute} palette {get_key_logscore(intervals)}; "
            "hide atoms;"
            "show cartoons;"
            f"key {get_key_logscore(intervals)} :{intervals[0]} :{intervals[1]} :{intervals[2]} :{intervals[3]} :{intervals[4]} pos 0.35,0.03 fontSize 4 size 0.3,0.02;"
            f"2dlabels create label text '{labels[attribute]}' size 6 color darkred xpos 0.35 ypos 0.065;"
            f"2dlabels create title text '{gene} - {uni_id}-F{f} ' size 6 color darkred xpos 0.35 ypos 0.93;"
            "lighting soft;"
            "graphics silhouettes true width 1.3;"
            "zoom;"
        )
    
    transparent_bg = " transparentBackground  true" if transparent_bg else ""
    
    if clusters is not None:
        if len(clusters) > 0:
            for pos in clusters:
                chimerax_command += f"marker #10 position :{pos} color #dacae961 radius 5.919;"
            chimerax_command += f"save {chimera_output_path}/{cohort}_{i}_{gene}_{attribute}_clusters.png pixelSize {pixelsize} supersample 3{transparent_bg};exit\""
    else:
        chimerax_command += f"save {chimera_output_path}/{cohort}_{i}_{gene}_{attribute}.png pixelSize {pixelsize} supersample 3{transparent_bg};exit\""
    
    return chimerax_command
            
            
def generate_chimerax_plot(output_dir,
                            gene_result_path,
                            pos_result_path,
                            datasets_dir,
                            seq_df_path,
                            cohort,
                            pixel_size,
                            cluster_ext,
                            fragmented_proteins,
                            palette,
                            transparent_bg,
                            chimerax_bin):

    chimera_out_path = f"{output_dir}/chimerax"
    chimera_attr_path = f"{chimera_out_path}/attributes"
    chimera_plots_path =f"{chimera_out_path}/plots"

    for path in [chimera_out_path, chimera_attr_path, chimera_plots_path]:
        if not os.path.exists(path):
            os.makedirs(path)

    seq_df = pd.read_csv(seq_df_path, sep="\t")
    result_genes = pd.read_csv(gene_result_path)
    result = pd.read_csv(pos_result_path)
    if "Ratio_obs_sim" in result.columns:
        result = result.rename(columns={"Ratio_obs_sim" : "Score_obs_sim"})
    result["Logscore_obs_sim"] = np.log(result["Score_obs_sim"])

    # Do this for each gene you want to plot (maybe top 30 or significant ones)
    genes = result_genes[result_genes["C_gene"] == 1].Gene.unique()
    if len(genes) > 0:
        for i, gene in enumerate(genes):
            logger.info(f"Processing {gene}")
            
            # Attribute files
            logger.info("Generating attribute files..")
            result_gene = result[result["Gene"] == gene]
            if cluster_ext:
                clusters = result_gene[result_gene.C == 1].Pos.values
            else:
                clusters = result_gene[(result_gene.C == 1) & (result_gene.C_ext == 0)].Pos.values            
            len_gene = len(seq_df[seq_df["Gene"] == gene].Seq.values[0])
            result_gene = pd.DataFrame({"Pos" : range(1, len_gene+1)}).merge(
                result_gene[["Pos", "Mut_in_res", "Mut_in_vol", "Score_obs_sim", "Logscore_obs_sim"]], on="Pos", how="left")

            logger.info("Generating 3D images..")
            uni_id, f = seq_df[seq_df["Gene"] == gene][["Uniprot_ID", "F"]].values[0]
            pdb_path = f"{datasets_dir}/pdb_structures/AF-{uni_id}-F{f}-model_v4.pdb"
            
            labels = {"mutres" : "Mutations in residue ", 
                    "mutvol" : "Mutations in volume ",
                    "score" : "   Clustering score ",
                    "logscore" : "log(Clustering score) "}
            cols = {"mutres" : "Mut_in_res", 
                    "mutvol" : "Mut_in_vol",
                    "score" : "Score_obs_sim",
                    "logscore" : "Logscore_obs_sim"}                                            
            
            if fragmented_proteins == False:
                if f != 1:
                    logger.info(f"Skipping {gene} ({uni_id}-F{f})..")
                    continue
                
            if os.path.exists(pdb_path):
                
                for attribute in ["mutres", "mutvol", "score", "logscore"]:
                    
                    attr_file_path = f"{chimera_attr_path}/{gene}_{attribute}.defattr"
                    create_attribute_file(path_to_file=attr_file_path,
                                        df=result_gene.dropna(),
                                        attribute_col=cols[attribute],
                                        attribute_name=attribute)
                    
                    attribute_vector = result_gene[cols[attribute]]
                    intervals = get_intervals(attribute_vector, attribute)
 
                    chimerax_command = get_chimerax_command(chimerax_bin, 
                                                            pdb_path, 
                                                            chimera_plots_path, 
                                                            attr_file_path, 
                                                            attribute, 
                                                            intervals, 
                                                            gene,
                                                            uni_id,
                                                            labels,
                                                            i,
                                                            f,
                                                            cohort,
                                                            pixelsize=pixel_size,
                                                            palette=palette,
                                                            transparent_bg=transparent_bg)
                    subprocess.run(chimerax_command, shell=True)
                    logger.debug(chimerax_command)
                    
                    if attribute == "score" or attribute == "logscore":
                        chimerax_command = get_chimerax_command(chimerax_bin, 
                                                                pdb_path, 
                                                                chimera_plots_path, 
                                                                attr_file_path, 
                                                                attribute, 
                                                                intervals, 
                                                                gene,
                                                                uni_id,
                                                                labels,
                                                                i,
                                                                f,
                                                                cohort,
                                                                clusters=clusters,
                                                                pixelsize=pixel_size,
                                                                palette=palette,
                                                                transparent_bg=transparent_bg)
                    subprocess.run(chimerax_command, shell=True)
                    logger.debug(chimerax_command)
            else:
                logger.info(f"PDB path missing: {pdb_path}")
    else:
        logger.info("Nothing to plot!")
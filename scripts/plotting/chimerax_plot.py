"""
Use ChimeraX to generate 3D structures colored by metrics
"""

import pandas as pd   
import numpy as np
import os 
import subprocess
import math
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


def get_palette(intervals, type="diverging"):
    
    # Diverging palette
    if type == "diverging":
        return f"{intervals[0]},#0571B0:{intervals[1]},#92C5DE:{intervals[2]},white:{intervals[3]},#F4A582:{intervals[4]},#CA0020"
    
    # Sequential palette
    else:
        return f"{intervals[0]},#FFFFB2:{intervals[1]},#FECC5C:{intervals[2]},#FD8D3C:{intervals[3]},#F03B20:{intervals[4]},#BD0026"


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
                         transparent_bg=False):
    
    palette = get_palette(intervals, type="diverging") if attribute == "logscore" else get_palette(intervals, type="sequential")
    transparent_bg = " transparentBackground  true" if transparent_bg else ""
    
    chimerax_command = (
        f"{chimerax_bin} --nogui --offscreen --silent --cmd "
        f"\"open {pdb_path}; "
        "set bgColor white; "
        "color lightgray; "
        f"open {attr_file_path}; "
        f"color byattribute {attribute} palette {palette}; "
        f"key {palette} :{intervals[0]} :{intervals[1]} :{intervals[2]} :{intervals[3]} :{intervals[4]} pos 0.35,0.03 fontSize 4 size 0.3,0.02;"
        f"2dlabels create label text '{labels[attribute]}' size 6 color darkred xpos 0.34 ypos 0.065;"
        f"2dlabels create title text '{gene} - {uni_id}-F{f} ' size 6 color darkred xpos 0.35 ypos 0.93;"
        "hide atoms;"
        "show cartoons;"
        "lighting soft;"
        "graphics silhouettes true width 1.3;"
        "zoom;"
    )
    
    if clusters is not None and len(clusters) > 0:
        for pos in clusters:
            chimerax_command += f"marker #10 position :{pos} color #dacae961 radius 5.919;"
        cluster_tag = "_clusters"
    else:
        cluster_tag = ""
    
    output_path = os.path.join(chimera_output_path, f"{cohort}_{i}_{gene}_{attribute}{cluster_tag}.png")
    chimerax_command += f"save {output_path} pixelSize {pixelsize} supersample 3{transparent_bg};exit\""
    
    return chimerax_command
            

def generate_chimerax_plot(output_dir,
                            gene_result_path,
                            pos_result_path,
                            datasets_dir,
                            seq_df_path,
                            cohort,
                            max_genes,
                            pixel_size,
                            cluster_ext,
                            fragmented_proteins,
                            transparent_bg,
                            chimerax_bin):

    seq_df = pd.read_csv(seq_df_path, sep="\t")
    gene_result = pd.read_csv(gene_result_path)
    result = pd.read_csv(pos_result_path)
    if "Ratio_obs_sim" in result.columns:
        result = result.rename(columns={"Ratio_obs_sim" : "Score_obs_sim"})
    result["Logscore_obs_sim"] = np.log(result["Score_obs_sim"])

    # Process each gene
    genes = gene_result[gene_result["C_gene"] == 1].Gene.unique()
    if len(genes) > 0:
        
        chimera_out_path = os.path.join(output_dir, f"{cohort}.chimerax")
        chimera_attr_path = os.path.join(chimera_out_path, "attributes")
        chimera_plots_path = os.path.join(chimera_out_path, "plots")
        for path in [chimera_out_path, chimera_attr_path, chimera_plots_path]:
            os.makedirs(path, exist_ok=True)
                
        for i, gene in enumerate(genes[:max_genes]):
            logger.info(f"Processing {gene}")
            
            # Attribute files
            logger.debug("Preprocessing for attribute files..")
            result_gene = result[result["Gene"] == gene]
            if cluster_ext:
                clusters = result_gene[result_gene.C == 1].Pos.values
            else:
                clusters = result_gene[(result_gene.C == 1) & (result_gene.C_ext == 0)].Pos.values            
            len_gene = len(seq_df[seq_df["Gene"] == gene].Seq.values[0])
            result_gene = pd.DataFrame({"Pos" : range(1, len_gene+1)}).merge(
                result_gene[["Pos", 
                             "Mut_in_res", 
                             "Mut_in_vol", 
                             "Score_obs_sim", 
                             "Logscore_obs_sim"]], on="Pos", how="left")

            uni_id, f = seq_df[seq_df["Gene"] == gene][["Uniprot_ID", "F"]].values[0]
            pdb_path = os.path.join(datasets_dir, "pdb_structures", f"AF-{uni_id}-F{f}-model_v4.pdb")
            
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
                    logger.debug(f"Fragmented protein processing {fragmented_proteins}: Skipping {gene} ({uni_id}-F{f})..")
                    continue
                
            if os.path.exists(pdb_path):
                
                for attribute in ["mutres", "mutvol", "score", "logscore"]:
                    
                    logger.debug("Generating attribute files..")
                    attr_file_path = f"{chimera_attr_path}/{gene}_{attribute}.defattr"
                    create_attribute_file(path_to_file=attr_file_path,
                                        df=result_gene.dropna(),
                                        attribute_col=cols[attribute],
                                        attribute_name=attribute)
                    
                    attribute_vector = result_gene[cols[attribute]]
                    intervals = get_intervals(attribute_vector, attribute)
 
                    logger.debug("Generating 3D images..")
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
                                                                transparent_bg=transparent_bg)
                        subprocess.run(chimerax_command, shell=True)
                        logger.debug(chimerax_command)
                        
            else:
                logger.warning(f"PDB path missing: {pdb_path}")
    else:
        logger.info("Nothing to plot!")
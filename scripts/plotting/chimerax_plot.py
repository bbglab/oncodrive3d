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

CHIMERA_BIN = "/usr/bin/chimerax"
    
def create_attribute_file(path_to_file, 
                                  df, 
                                  attribute_col,
                                  pos_col="Pos",
                                  attribute_name="local_attribute",
                                  gene="Gene_name", 
                                  protein="Protein_name"):

    print(f"Saving {path_to_file}")
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


def get_intervals(attribute_vector, attribute, n=5):
    
    max_value = int(attribute_vector.max())
    intervals = [i*int(max_value / 4) for i in range(5)]
    intervals[0] = round_to_first_nonzero(attribute_vector.min()) if attribute == "score" else 1
    intervals[4] = max_value

    if attribute == "logscore":
        min_value = np.round(np.min(attribute_vector), 1)
        max_value = np.round(np.max(attribute_vector), 1)
        min_value = round(min_value) if is_float_an_integer(min_value) else min_value
        max_value = round(max_value) if is_float_an_integer(max_value) else max_value   
        
        if max_value <= 0:
            pos_val = 0
        else:
            pos_val = round(max_value / 2)
            
        if min_value >= 0:
            neg_val = 0
        else:
            neg_val = round(min_value / 2)
        
        intervals = [min_value, neg_val, 0, pos_val, max_value]
        
    return intervals


def get_key_logscore(intervals):
    return f"{intervals[0]},#0571B0:{intervals[1]},#92C5DE:{intervals[2]},white:{intervals[3]},#F4A582:{intervals[4]},#CA0020"


def get_chimerax_command(chimera_bin, 
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
        f"{chimera_bin} --nogui --offscreen --cmd "
        f"\"open {pdb_path}; "
        "set bgColor white; "
        "color lightgray; "
        f"open {attr_file_path}; "
        f"color byattribute {attribute} palette {palette}; "
        "hide atoms;"
        "show cartoons;"
        f"key {palette} :{intervals[0]} :{intervals[1]}  :{intervals[2]} :{intervals[3]} :{intervals[4]} pos 0.35,0.03 fontSize 4 size 0.3,0.02;"
        f"2dlabels create label text '{labels[attribute]}' size 6 color darkred xpos 0.365 ypos 0.065;"
        f"2dlabels create title text '{gene} - {uni_id}-F{f} ' size 6 color darkred xpos 0.375 ypos 0.93;"
        "lighting soft;"
        "graphics silhouettes true width 1.3;"
        "zoom;"
    )
    
    if attribute == "logscore":
        chimerax_command = (
            f"{chimera_bin} --nogui --offscreen --cmd "
            f"\"open {pdb_path}; "
            "set bgColor white; "
            "color lightgray; "
            f"open {attr_file_path}; "
            f"color byattribute {attribute} palette {get_key_logscore(intervals)}; "
            "hide atoms;"
            "show cartoons;"
            f"key {get_key_logscore(intervals)} :{intervals[0]} :{intervals[1]} :{intervals[2]} :{intervals[3]} :{intervals[4]} pos 0.35,0.03 fontSize 4 size 0.3,0.02;"
            f"2dlabels create label text '{labels[attribute]}' size 6 color darkred xpos 0.365 ypos 0.065;"
            f"2dlabels create title text '{gene} - {uni_id}-F{f} ' size 6 color darkred xpos 0.375 ypos 0.93;"
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
            

@click.command(context_settings=dict(help_option_names=['-h', '--help']),
               help="Generate 3D plots using CHimeraX.")
@click.option("-o", "--output_dir", 
              help="Directory where to save the plots", type=str, required=True)
@click.option("-g", "--gene_result_path", 
              help="Path to genes-level O3D result", type=click.Path(exists=True), required=True)
@click.option("-p", "--pos_result_path", 
              help="Path to positions-level O3D result", type=click.Path(exists=True), required=True)
@click.option("-d", "--datasets_dir", 
              help="Path to datasets", type=click.Path(exists=True), required=True)
@click.option("-s", "--seq_df_path", 
              help="Path to sequences dataframe", type=click.Path(exists=True), required=True)
@click.option("-c", "--cohort", 
              help="Cohort name", default="")
@click.option("--pixel_size", help="Pixel size (smaller value is larger number of pixels)", type=float, default=0.1)
@click.option("--cluster_ext", help="Include extended clusters", is_flag=True)
@click.option("--fragmented_proteins", help="Include fragmented proteins", is_flag=True)
@click.option("--palette", help="Color palette", type=str, default="YlOrRd-5")
@click.option("--transparent_bg", help="Set background as transparent", type=str, is_flag=True)
def main(output_dir,
         gene_result_path,
         pos_result_path,
         datasets_dir,
         seq_df_path,
         cohort,
         pixel_size,
         cluster_ext,
         fragmented_proteins,
         palette,
         transparent_bg):

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
    for i, gene in enumerate(genes):
        print("\n> Processing", gene)
        
        # Attribute files
        print("Generating attribute files..")
        result_gene = result[result["Gene"] == gene]
        if cluster_ext:
            clusters = result_gene[result_gene.C == 1].Pos.values
        else:
            clusters = result_gene[(result_gene.C == 1) & (result_gene.C_ext == 0)].Pos.values            
        len_gene = len(seq_df[seq_df["Gene"] == gene].Seq.values[0])
        result_gene = pd.DataFrame({"Pos" : range(1, len_gene+1)}).merge(
            result_gene[["Pos", "Mut_in_res", "Mut_in_vol", "Score_obs_sim", "Logscore_obs_sim"]], on="Pos", how="left")

        print("Generating 3D images..")
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
                print(f"Skipping {gene} ({uni_id}-F{f})..")
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
                
                chimerax_command = get_chimerax_command(CHIMERA_BIN, 
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
                
                if attribute == "score" or attribute == "logscore":
                    chimerax_command = get_chimerax_command(CHIMERA_BIN, 
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
        else:
            print(f"PDB path missing: {pdb_path}")
            
            
if __name__ == "__main__":
    main()
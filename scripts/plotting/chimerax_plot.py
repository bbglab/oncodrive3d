"""
Use ChimeraX to generate 3D structures colored by metrics
"""

import os
import shlex
import subprocess
import math

import pandas as pd   
import numpy as np
import daiquiri

from scripts import __logger_name__
from scripts.plotting.utils import cap_inf_scores, detect_af_version

logger = daiquiri.getLogger(__logger_name__ + ".plotting.chimerax_plot")


# Paths are interpolated unquoted into the --cmd script (ChimeraX rejects quoted
# paths in a ';'-separated --cmd on some versions), so reject characters that
# would split a path (space) or inject extra commands (quote, newline).
_CHIMERAX_UNSAFE_PATH_CHARS = ('"', ' ', '\n', '\r')


def _validate_chimerax_path(path, role):
    """Reject paths with characters that would break the ChimeraX --cmd script.

    Real-world paths here are normally clean; this just surfaces a clear failure
    instead of a silently broken or injected command.
    """
    text = str(path)
    bad = [c for c in _CHIMERAX_UNSAFE_PATH_CHARS if c in text]
    if bad:
        raise ValueError(
            f"{role} path contains characters unsafe for the ChimeraX command "
            f"script ({', '.join(repr(c) for c in bad)}): {text!r}"
        )


def _run_chimerax(command, label):
    """Run a ChimeraX subprocess. Failures become warnings so one bad image
    doesn't abort the loop while the user still sees something went wrong.

    Paths inside the --cmd script are unquoted, so `get_chimerax_command`
    validates them to be free of spaces/quotes/newlines that would break the
    ChimeraX parser.
    """
    # Discard stdout (never read) but keep stderr for the failure warning;
    # capture_output=True would buffer both into memory across many genes.
    try:
        result = subprocess.run(
            command,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
        )
    except OSError as exc:
        # E.g. FileNotFoundError when chimerax_bin is missing or not executable.
        # Treat it like a failed run so the loop continues across other genes.
        logger.warning(f"ChimeraX could not be executed for {label}: {exc}")
        return None
    if result.returncode != 0:
        logger.warning(
            f"ChimeraX failed for {label} "
            f"(rc={result.returncode}): {result.stderr.strip()}"
        )
    return result


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
                         colored_positions=None,
                         pixelsize=0.1,
                         transparent_bg=False):

    palette = get_palette(intervals, type="diverging") if attribute == "logscore" else get_palette(intervals, type="sequential")
    transparent_bg_suffix = " transparentBackground  true" if transparent_bg else ""

    # Validate before interpolating: paths go unquoted into the --cmd script.
    _validate_chimerax_path(pdb_path, "PDB")
    _validate_chimerax_path(attr_file_path, "attribute file")

    chimerax_script = (
        f"open {pdb_path}; "
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

    # Render mutated (= colored) residues as spheres so their colour is actually
    # visible on the cartoon. Without this, the colour mapping is washed out by
    # the cartoon ribbon. ~sel afterwards clears the green selection markers
    # from the saved image.
    if colored_positions is not None and len(colored_positions) > 0:
        pos_selector = ",".join(str(int(p)) for p in sorted(set(colored_positions)))
        chimerax_script += (
            f"sel :{pos_selector};"
            "style sel sphere;"
            "hide sel cartoons;"
            "show sel atoms;"
            "~sel;"
        )

    if clusters is not None and len(clusters) > 0:
        for pos in clusters:
            chimerax_script += f"marker #10 position :{pos} color #dacae961 radius 5.919;"
        cluster_tag = "_clusters"
    else:
        cluster_tag = ""

    output_path = os.path.join(chimera_output_path, f"{cohort}_{i}_{gene}_{attribute}{cluster_tag}.png")
    _validate_chimerax_path(output_path, "output")
    chimerax_script += f'save {output_path} pixelSize {pixelsize} supersample 3{transparent_bg_suffix};exit'

    # argv list for subprocess.run(shell=False) — the OS shell never sees it.
    return [chimerax_bin, "--nogui", "--offscreen", "--silent", "--cmd", chimerax_script]
            

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
                            chimerax_bin,
                            af_version):

    seq_df = pd.read_csv(seq_df_path, sep="\t")
    gene_result = pd.read_csv(gene_result_path)
    result = pd.read_csv(pos_result_path)
    if result.empty:
        logger.warning("Empty position-level result; nothing to plot.")
        return
    if "Ratio_obs_sim" in result.columns:
        result = result.rename(columns={"Ratio_obs_sim" : "Score_obs_sim"})
    result = cap_inf_scores(result)
    result["Logscore_obs_sim"] = np.log(result["Score_obs_sim"])

    # Detect the AlphaFold version actually present in the dataset; the
    # user-supplied --af_version is used as a tiebreaker only.
    af_version = detect_af_version(datasets_dir, requested_version=af_version)

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
            pdb_path_base = os.path.join(datasets_dir, "pdb_structures", f"AF-{uni_id}-F{f}-model_v{af_version}.pdb")
            pdb_candidates = [pdb_path_base, f"{pdb_path_base}.gz"]
            pdb_path = next((candidate for candidate in pdb_candidates if os.path.exists(candidate)), None)
            
            labels = {"mutres" : "Mutations in residue ", 
                      "mutvol" : "Mutations in volume ",
                      "score" : "   Clustering score ",
                      "logscore" : "log(Clustering score) "}
            cols = {"mutres" : "Mut_in_res", 
                    "mutvol" : "Mut_in_vol",
                    "score" : "Score_obs_sim",
                    "logscore" : "Logscore_obs_sim"}                                            
            
            if fragmented_proteins == False:
                if str(f) != "1":
                    logger.debug(f"Fragmented protein processing {fragmented_proteins}: Skipping {gene} ({uni_id}-F{f})..")
                    continue
                
            if pdb_path:

                # Mutated rows: the residues written to the .defattr file and
                # thus coloured. Reused below so the spheres match the colours.
                result_gene_mutated = result_gene.dropna()
                colored_positions = result_gene_mutated["Pos"].astype(int).tolist()

                for attribute in ["mutres", "mutvol", "score", "logscore"]:

                    logger.debug("Generating attribute files..")
                    attr_file_path = f"{chimera_attr_path}/{gene}_{attribute}.defattr"
                    create_attribute_file(path_to_file=attr_file_path,
                                        df=result_gene_mutated,
                                        attribute_col=cols[attribute],
                                        attribute_name=attribute)

                    attribute_vector = result_gene[cols[attribute]]
                    intervals = get_intervals(attribute_vector, attribute)

                    logger.debug("Generating 3D images..")
                    try:
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
                                                                colored_positions=colored_positions,
                                                                pixelsize=pixel_size,
                                                                transparent_bg=transparent_bg)
                    except ValueError as exc:
                        logger.warning(f"Skipping {gene} ({attribute}): {exc}")
                        continue
                    _run_chimerax(chimerax_command, f"{gene} ({attribute})")
                    logger.debug(shlex.join(chimerax_command))

                    if attribute == "score" or attribute == "logscore":
                        try:
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
                                                                    colored_positions=colored_positions,
                                                                    pixelsize=pixel_size,
                                                                    transparent_bg=transparent_bg)
                        except ValueError as exc:
                            logger.warning(f"Skipping {gene} ({attribute}, clusters): {exc}")
                            continue
                        _run_chimerax(chimerax_command, f"{gene} ({attribute}, clusters)")
                        logger.debug(shlex.join(chimerax_command))
                        
            else:
                tried_files = ', '.join(os.path.basename(path) for path in pdb_candidates)
                logger.warning(f"Structure not found for {uni_id}-F{f} (AlphaFold v{af_version}). Tried: {tried_files}")
    else:
        logger.info("Nothing to plot!")

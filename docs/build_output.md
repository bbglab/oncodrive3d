# Building Datasets Output

After `oncodrive3d build-datasets` completes (or after extracting the pre-built datasets [from Zenodo](https://zenodo.org/records/21031511)), the build folder contains:

```text
<build_folder>/
├── pdb_structures/        # AlphaFold PDB files (downloaded bundle + any custom-supplied structures)
├── pae/                   # AlphaFold Predicted Aligned Error files (.npy, one per UniProt protein — fragmented proteins use the F1 fragment's PAE); absent or partial when PAE is unavailable for the chosen AF version
├── prob_cmaps/            # Contact probability maps (.npy, one per protein/fragment); F2+ fragments and proteins without PAE fall back to binary contact maps
├── confidence.tsv         # Per-residue pLDDT scores for every protein in the proteome
├── seq_for_mut_prob.tsv   # Per-protein metadata: HUGO symbol, HGNC/Ensembl/UniProt IDs, protein and DNA sequence, exon coordinates, trinucleotide contexts
├── biomart_metadata.tsv   # Ensembl BioMart canonical transcript metadata; used to prioritize when multiple transcripts map to one structure
└── log/                   # Logs from the build run
```

With `--mane` or `--mane_only`, the folder also contains assets from the AlphaFold MANE overlap bundle (prefixed with `mane_`):

```text
├── mane_summary.txt.gz                # NCBI MANE release summary (cached after first download)
├── mane_refseq_prot_to_alphafold.csv  # RefSeq → AlphaFold MANE bundle mapping (consumed by the MANE preprocessing toolkit)
├── mane_missing.csv                   # MANE entries lacking an AlphaFold MANE structure (consumed by the MANE preprocessing toolkit)
└── mane_readme.txt                    # README from the AlphaFold MANE overlap bundle
```

> [!NOTE]
> If `--custom_pae_dir` is provided, `pae/` is replaced with the contents of that directory. When PAE files are unavailable (e.g., AF DB v4 after 2025) and no custom PAE is supplied, the build proceeds with binary contact maps in `prob_cmaps/`.

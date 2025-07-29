# Output

After successfully completing this step, the build folder must include the following files and directories:

- **pae/**: Directory including the AlphaFold predicted aligned error (PAE) for any protein of the proteome with a length lower than 2700 amino acids

- **prob_cmaps/**: Directory including the contact probability map (pCMAPs) for any protein of the proteome

- **confidence.csv**: CSV file including per-residue predicted local distance difference test (pLDDT) score for any protein of the proteome

- **seq_for_mut_prob.csv**: CSV file including HUGO symbol, Uniprot ID, DNA and protein sequences for any proteine of the proteome

- **pdb_structures/**: Directory containing all PDB structure files—both those downloaded from the AlphaFold database and any custom in-house predicted structures (if provided).

- **log/**: Directory containing log files produced during the build-datasets execution.

- **biomart_metadata.csv**: Metadata file from Ensembl BioMart used to prioritize canonical transcripts when multiple transcripts map to the same protein structure.

...
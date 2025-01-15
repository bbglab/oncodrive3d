# Output

After successfully completing this step, the build folder must include the following files and directories:

- **pae/**: directory including the AlphaFold predicted aligned error (PAE) for any protein of the proteome with a length lower than 2700 amino acids

- **prob_cmaps/**: directory including the contact probability map (pCMAPs) for any protein of the proteome

- **confidence.csv**: CSV file including per-residue predicted local distance difference test (pLDDT) score for any protein of the proteome

- **seq_for_mut_prob.csv**: CSV file including HUGO symbol, Uniprot ID, DNA and protein sequences for any proteine of the proteome

...
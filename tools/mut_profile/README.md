# Mutational profile helpers

## `get_regions_file.py`

Generates the genomic **regions file** (`CHROMOSOME`, `START`, `END` for the
autosomes plus `X`/`Y`) needed to compute a cohort mutational profile.

```bash
python get_regions_file.py <genome> <kmer> > regions.tsv
```

- `<genome>`: one of `hg18`, `hg19`, `hg38`, `mm10`, `mm39`.
- `<kmer>`: context width (e.g. `3` for trinucleotide). Each chromosome is
  trimmed by `kmer // 2` at both ends so every position has full context.

Requires [`bgreference`](https://pypi.org/project/bgreference/).

See [Create the Mutation Profile with BGSignature](../../docs/run_input_output.md#create-the-mutation-profile-with-bgsignature)
for how this file feeds into building the profile passed to `oncodrive3d run -p`.

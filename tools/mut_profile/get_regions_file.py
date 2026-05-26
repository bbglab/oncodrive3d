import sys

from bgreference import refseq


GENOME_AUTOSOMES = {
    "hg18": 22, "hg19": 22, "hg38": 22,
    "mm10": 19, "mm39": 19,
}


def compute_sizes(genome, kmer):
    if genome not in GENOME_AUTOSOMES:
        raise ValueError(
            f"Unsupported genome: {genome!r}. Supported: {sorted(GENOME_AUTOSOMES)}"
        )
    chroms = [str(i) for i in range(1, GENOME_AUTOSOMES[genome] + 1)] + ['X', 'Y']
    sizes = []
    for chr_ in chroms:
        seq = refseq(genome, chr_, start=(1 + kmer // 2), size=None)
        sizes.append(tuple(map(str, (chr_, 1 + kmer // 2, len(seq) - kmer // 2))))
    return sizes


def write(sizes):
    print('\t'.join(('CHROMOSOME', 'START', 'END')))
    for s in sizes:
        print('\t'.join(s))


if __name__ == '__main__':
    genome = sys.argv[1]
    kmer = int(sys.argv[2])
    write(compute_sizes(genome, kmer))

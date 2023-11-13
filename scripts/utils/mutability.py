"""
This module contains the methods associated with the
mutabilities that are assigned to the mutations.

The mutabilities are read from a file.
The file must be compressed using bgzip, and then indexed using tabix.
$ bgzip ..../all_samples.mutability_per_site.tsv
$ tabix -b 2 -e 2 ..../all_samples.mutability_per_site.tsv.gz

"""
import logging
import gzip
import mmap
import json
import os
import struct
from collections import defaultdict, namedtuple
from typing import List

import bgdata
import numpy as np
import tabix

# from scripts import __logger_name__
# from oncodrivefml.reference import get_ref_triplet

# logger = logging.getLogger(__logger_name__)

transcribe = {"A":"T", "C":"G", "G":"C", "T":"A"}

MutabilityValue = namedtuple('MutabilityValue', ['ref', 'alt', 'value'])
"""
Tuple that contains the reference, the alteration, the mutability value

Parameters:
    ref (str): reference base
    alt (str): altered base
    value (float): mutability value of that substitution
"""

mutabilities_reader = None


class ReaderError(Exception):

    def __init__(self, msg):
        self.message = msg


class ReaderGetError(ReaderError):
    def __init__(self, chr, start, end):
        self.message = 'Error reading chr: {} start: {} end: {}'.format(chr, start, end)

class MutabilityTabixReader:

    def __init__(self, conf):
        self.file = conf['file']
        self.conf_chr_prefix = conf['chr_prefix']
        self.ref_pos = conf['ref']
        self.alt_pos = conf['alt']
        self.pos_pos = conf['pos']
        self.mutability_pos = conf['mutab']
        self.element_pos = None

    def __enter__(self):
        self.tb = tabix.open(self.file)
        self.index_errors = 0
        self.elements_errors = 0
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.index_errors > 0 or self.elements_errors > 0:
            raise ReaderError('{} index errors and {} discrepancies between the expected and retrieved element'.format(self.index_errors, self.elements_errors))
        return True

    def _read_row(self, row):
        mutability = float(row[self.mutability_pos])
        ref = None if self.ref_pos is None else row[self.ref_pos]
        alt = None if self.alt_pos is None else row[self.alt_pos]
        pos = None if self.pos_pos is None else int(row[self.pos_pos])
        element = None if self.element_pos is None else row[self.element_pos]
        return (mutability, ref, alt, pos), element

    def get(self, chromosome, start, stop, element=None):
        try:
            for row in self.tb.query("{}{}".format(self.conf_chr_prefix, chromosome), start, stop):
                try:
                    r = self._read_row(row)
                except IndexError:
                    self.index_errors += 1
                    continue
                else:
                    if self.element_pos is not None and element is not None and r[1] != element:
                        self.elements_errors += 1
                        continue
                    yield r[0]
        except tabix.TabixError:
            raise ReaderGetError(chromosome, start, stop)


def init_mutabilities_module(conf):
    global mutabilities_reader
    # TODO add an else case or fix this function
    mutabilities_reader = MutabilityTabixReader(conf)


class Mutabilities(object):
    """

    Args:
        element (str): element ID
        segments (list): list of the segments associated to the element
        config (dict): configuration

    Attributes:
        mutabilities_by_pos (dict): for each positions get all possible changes

            .. code-block:: python

                    { position:
                        [
                            MutabilityValue(
                                ref,
                                alt_1,
                                value
                            ),
                            MutabilityValue(
                                ref,
                                alt_2,
                                value
                            ),
                            MutabilityValue(
                                ref,
                                alt_3,
                                value
                            )
                        ]
                    }
    """

    def __init__(self, element: str, chromosome:str, segments: list, gene_len: int, gene_strand: bool, config: dict):

        
        self.element = element
        self.chromosome = chromosome
        self.segments = segments
        self.gene_length = gene_len
        self.gene_strand = gene_strand
        # [{'CHROMOSOME': '10', 'START': 87864532, 'END': 87864548, 'ELEMENT': 'PTEN', 'SEGMENT': 'exon01', 'LINE': 77},
        #     {'CHROMOSOME': '10', 'START': 87894024, 'END': 87894109, 'ELEMENT': 'PTEN', 'SEGMENT': 'exon02', 'LINE': 78},
        #     {'CHROMOSOME': '10', 'START': 87925512, 'END': 87925557, 'ELEMENT': 'PTEN', 'SEGMENT': 'exon03', 'LINE': 79},
        #     {'CHROMOSOME': '10', 'START': 87931045, 'END': 87931089, 'ELEMENT': 'PTEN', 'SEGMENT': 'exon04', 'LINE': 80},
        #     {'CHROMOSOME': '10', 'START': 87933012, 'END': 87933022, 'ELEMENT': 'PTEN', 'SEGMENT': 'exon05', 'LINE': 81},
        #     {'CHROMOSOME': '10', 'START': 87933223, 'END': 87933251, 'ELEMENT': 'PTEN', 'SEGMENT': 'exon05', 'LINE': 82},
        #     {'CHROMOSOME': '10', 'START': 87952117, 'END': 87952135, 'ELEMENT': 'PTEN', 'SEGMENT': 'exon06', 'LINE': 83},
        #     {'CHROMOSOME': '10', 'START': 87952220, 'END': 87952259, 'ELEMENT': 'PTEN', 'SEGMENT': 'exon06', 'LINE': 84},
        #     {'CHROMOSOME': '10', 'START': 87957852, 'END': 87957890, 'ELEMENT': 'PTEN', 'SEGMENT': 'exon07', 'LINE': 85},
        #     {'CHROMOSOME': '10', 'START': 87957961, 'END': 87958019, 'ELEMENT': 'PTEN', 'SEGMENT': 'exon07', 'LINE': 86},
        #     {'CHROMOSOME': '10', 'START': 87960893, 'END': 87960907, 'ELEMENT': 'PTEN', 'SEGMENT': 'exon08', 'LINE': 87},
        #     {'CHROMOSOME': '10', 'START': 87965286, 'END': 87965296, 'ELEMENT': 'PTEN', 'SEGMENT': 'exon09', 'LINE': 88}
        # ]
        
        # mutability configuration
        self.conf_file = config['file']
        # self.conf_mutability = config['mutab']
        self.conf_chr = config['chr']
        self.conf_chr_prefix = config['chr_prefix']
        self.conf_ref = config['ref']
        self.conf_alt = config['alt']
        self.conf_pos = config['pos']
        self.conf_element = config['element'] if 'element' in config.keys() else None
        self.conf_extra = config['extra'] if 'extra' in config.keys() else None

        # mutabilities to load
        self.mutabilities_by_pos = defaultdict(dict)


        # Initialize background mutabilities
        self._load_mutabilities()

    def get_mutability_by_position(self, position: int):
        """
        Get all MutabilityValue objects that are associated with that position

        Args:
            position (int): position

        Returns:
            :obj:`list` of :obj:`MutabilityValue`: list of all MutabilityValue related to that position

        """
        return self.mutabilities_by_pos.get(position, [])

    def get_all_positions(self) -> List[int]:
        """
        Get all positions in the element

        Returns:
            :obj:`list` of :obj:`int`: list of positions

        """
        return self.mutabilities_by_pos.keys()

    def _load_mutabilities(self):
        """
        For each position get all possible substitutions and for each
        obtains the assigned mutability

        Returns:
            dict: for each positions get a list of MutabilityValue
            (see :attr:`mutabilities_by_pos`)
        """
        segment_starting_pos = 0
        start = 0 if self.gene_strand else 1
        end = 1 if self.gene_strand else 0
        update_pos = 1 if self.gene_strand else -1
        prev_pos = None
        try:
            with mutabilities_reader as reader:
                for region in self.segments:
                    # region  # this would be an exon
                    try:
                        segment_len = region[end] - region[start] + 1
                        cdna_pos = segment_starting_pos if self.gene_strand else segment_starting_pos + segment_len
                        for row in reader.get(self.chromosome, region[start], region[end], self.element):
                            # every row is a site

                            mutability, ref, alt, pos = row
                            #print(mutability, ref, alt, pos, sep = "\t")

                            # ref_triplet = get_ref_triplet(region['CHROMOSOME'], pos - 1)
                            # ref = ref_triplet[1] if ref is None else ref
                            # if ref_triplet[1] != ref:
                            #     logger.warning("Background mismatch at position %d at '%s'", pos, self.element)


                            # if the current position is different from the previous
                            # update the cdna position accordingly to the strand
                            # and also update the value of prev_pos
                            if pos != prev_pos:
                                cdna_pos += update_pos
                                prev_pos = pos
                            
                            # since at protein level we are looking at the nucleotide 
                            # changes of the translated codons we store them as they will be queried later
                            if not self.gene_strand:
                                alt = transcribe[alt]

                            # add the mutability
                            self.mutabilities_by_pos[cdna_pos][alt] = mutability

                        segment_starting_pos = cdna_pos if self.gene_strand else cdna_pos + segment_len

                    except ReaderError as e:
                        logger.warning(e.message)
                        continue
        except ReaderError as e:
            logger.warning("Reader error: %s. Regions being analysed %s", e.message, self.segments)



if __name__ == "__main__":
    mutab_config = json.load(open('/home/fcalvet/Documents/dev/clustering_3d/test/normal_tests/mutability_config.json'))
    init_mutabilities_module(mutab_config)
    chrom = 17
    exons = eval("[(7676594, 7676521)]")
    seq_len = len("ATGGAGGAGCCGCAGTCAGATCCTAGCGTCGAGCCCCCTCTGAGTCAGGAAACATTTTCAGACCTATGGAAACT")
    # exons = eval("[(7676594, 7676521), (7676403, 7676382), (7676272, 7675994), (7675236, 7675053), (7674971, 7674859), (7674290, 7674181), (7673837, 7673701), (7673608, 7673535), (7670715, 7670609), (7669690, 7669612)]")
    tot_s_ex = 0
    for s, e in exons:
        tot_s_ex += np.sqrt((e-s)**2) + 1
    #print(tot_s_ex)
    # seq_len = len("ATGGAGGAGCCCCAGAGCGACCCCAGCGTGGAGCCCCCCCTGAGCCAGGAGACCTTCAGCGACCTGTGGAAGCTGCTGCCCGAGAACAACGTGCTGAGCCCCCTGCCCAGCCAGGCCATGGACGACCTGATGCTGAGCCCCGACGACATCGAGCAGTGGTTCACCGAGGACCCCGGCCCCGACGAGGCCCCCAGGATGCCCGAGGCCGCCCCCCCCGTGGCCCCCGCCCCCGCCGCCCCCACCCCCGCCGCCCCCGCCCCCGCCCCCAGCTGGCCCCTGAGCAGCAGCGTGCCCAGCCAGAAGACCTACCAGGGCAGCTACGGCTTCAGGCTGGGCTTCCTGCACAGCGGCACCGCCAAGAGCGTGACCTGCACCTACAGCCCCGCCCTGAACAAGATGTTCTGCCAGCTGGCCAAGACCTGCCCCGTGCAGCTGTGGGTGGACAGCACCCCCCCCCCCGGCACCAGGGTGAGGGCCATGGCCATCTACAAGCAGAGCCAGCACATGACCGAGGTGGTGAGGAGGTGCCCCCACCACGAGAGGTGCAGCGACAGCGACGGCCTGGCCCCCCCCCAGCACCTGATCAGGGTGGAGGGCAACCTGAGGGTGGAGTACCTGGACGACAGGAACACCTTCAGGCACAGCGTGGTGGTGCCCTACGAGCCCCCCGAGGTGGGCAGCGACTGCACCACCATCCACTACAACTACATGTGCAACAGCAGCTGCATGGGCGGCATGAACAGGAGGCCCATCCTGACCATCATCACCCTGGAGGACAGCAGCGGCAACCTGCTGGGCAGGAACAGCTTCGAGGTGAGGGTGTGCGCCTGCCCCGGCAGGGACAGGAGGACCGAGGAGGAGAACCTGAGGAAGAAGGGCGAGCCCCACCACGAGCTGCCCCCCGGCAGCACCAAGAGGGCCCTGCCCAACAACACCAGCAGCAGCCCCCAGCCCAAGAAGAAGCCCCTGGACGGCGAGTACTTCACCCTGCAGATCAGGGGCAGGGAGAGGTTCGAGATGTTCAGGGAGCTGAACGAGGCCCTGGAGCTGAAGGACGCCCAGGCCGGCAAGGAGCCCGGCGGCAGCAGGGCCCACAGCAGCCACCTGAAGAGCAAGAAGGGCCAGAGCACCAGCAGGCACAAGAAGCTGATGTTCAAGACCGAGGGCCCCGACAGCGAC")
    mutability_obj = Mutabilities("TP53", chrom, exons, seq_len, False, mutab_config)
    
    # for s, e in exons:
    #     for ii in range(min(s, e), max(s, e)+1):
    #         print(ii)

#    for key in mutability_obj.mutabilities_by_pos.keys():
#        print(key)
#        if len(mutability_obj.mutabilities_by_pos[key]) != 3:
#            print(mutability_obj.mutabilities_by_pos[key])

    for key in sorted(mutability_obj.mutabilities_by_pos):
        # print(key)
        print(key, mutability_obj.mutabilities_by_pos[key])
        # if len(mutability_obj.mutabilities_by_pos[key]) != 3:
        #     print(mutability_obj.mutabilities_by_pos[key])


    print(len(mutability_obj.mutabilities_by_pos))
    print(seq_len)
    mutability_dict = mutability_obj.mutabilities_by_pos
    
    # TODO raise an error here
    if len(mutability_dict) != seq_len:
        print("error")


    # TODO: reverse the nucleotides when working with the reverse strand
    
    
    # see which positions are lost in the process

    # 7676152 7676152
    # 7676153 7676153
    # 7676154 7676155
    # 7676155 7676156
    # 7676156 7676157
    # 7676157 7676158

    # 17      7676153 G       A       798.2016248293505
    # 17      7676153 G       C       128.15752499418232
    # 17      7676153 G       T       204.81320237872495
    # 17      7676155 G       A       233.83683455274522
    # 17      7676155 G       C       146.25715509711932
    # 17      7676155 G       T       162.44615179575143
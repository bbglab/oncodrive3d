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

from scripts import __logger_name__
from oncodrivefml.reference import get_ref_triplet

logger = logging.getLogger(__logger_name__)


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
        self.element_pos = conf['element']

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

    def __init__(self, element: str, segments: list, config: dict):

        self.element = element
        self.segments = segments
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
        self.conf_element = config['element']
        self.conf_extra = config['extra']

        # mutabilities to load
        self.mutabilities_by_pos = defaultdict(list)


        # Initialize background mutabilities
        self._load_mutabilities()

    def get_mutability_by_position(self, position: int) -> List[MutabilityValue]:
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

        try:
            with mutabilities_reader as reader:
                for region in self.segments:
                    try:
                        # region  # this would be an exon
                        for row in reader.get(region['CHROMOSOME'], region['START'], region['END'], self.element):
                            mutability, ref, alt, pos = row
                            ref_triplet = get_ref_triplet(region['CHROMOSOME'], pos - 1)
                            ref = ref_triplet[1] if ref is None else ref

                            if ref_triplet[1] != ref:
                                logger.warning("Background mismatch at position %d at '%s'", pos, self.element)

                            alts = alt if alt is not None and alt != '.' else 'ACGT'.replace(ref, '')

                            for a in alts:
                                self.mutabilities_by_pos[pos].append(MutabilityValue(ref, a, mutability))

                    except ReaderError as e:
                        logger.warning(e.message)
                        continue
        except ReaderError as e:
            logger.warning("Reader error: %s. Regions being analysed %s", e.message, self.segments)
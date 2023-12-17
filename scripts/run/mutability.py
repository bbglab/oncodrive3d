"""
This module contains the methods associated with the
mutabilities that are assigned to the mutations.

The mutabilities are read from a file.
The file must be compressed using bgzip, and then indexed using tabix.
$ bgzip ..../all_samples.mutability_per_site.tsv
$ tabix -b 2 -e 2 ..../all_samples.mutability_per_site.tsv.gz
"""

import logging
from collections import defaultdict, namedtuple
from typing import List

import tabix

# TODO uncomment these lines to make sure that the logger and everything is working

from scripts import __logger_name__
# from oncodrivefml.reference import get_ref_triplet

logger = logging.getLogger(__logger_name__ + ".run.mutability")

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

    def __init__(self, element: str, chromosome:str, segments: list, gene_len: int, gene_reverse_strand: bool, config: dict):

        
        self.element = element
        self.chromosome = chromosome
        self.segments = segments
        self.gene_length = gene_len
        self.reverse = True if float(gene_reverse_strand) == 1.0 else False

        # mutability configuration
        self.conf_file = config['file']
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
        cdna_pos = 0 
        starting_cdna_pos = 0
        start = 0 if not self.reverse else 1
        end = 1 if not self.reverse else 0
        update_pos = 1 if not self.reverse else -1
        try:
            with mutabilities_reader as reader:
                for region in self.segments:
                    # each region corresponds to an exon
                    try:
                        segment_len = region[end] - region[start] + 1
                        cdna_pos = starting_cdna_pos if not self.reverse else starting_cdna_pos + segment_len
                        starting_cdna_pos = int(cdna_pos)
                        prev_pos = region[start] - 1
                        # print(self.chromosome, region[start], region[end], self.element, segment_len, cdna_pos, prev_pos, starting_cdna_pos, update_pos)
                        for row in reader.get(self.chromosome, region[start], region[end], self.element):
                            # every row is a site

                            mutability, ref, alt, pos = row
                            # print(mutability, ref, alt, pos, sep = "\t")

                            # TODO decide if we want to do this check or not,
                            # I would say it is fine as it is now without the check

                            # ref_triplet = get_ref_triplet(region['CHROMOSOME'], pos - 1)
                            # ref = ref_triplet[1] if ref is None else ref
                            # if ref_triplet[1] != ref:
                            #     logger.warning("Background mismatch at position %d at '%s'", pos, self.element)


                            # if the current position is different from the previous
                            # update the cdna position accordingly to the strand
                            # and also update the value of prev_pos
                            if pos != prev_pos:
                                cdna_pos += update_pos
                                # print("changing position", pos, prev_pos, cdna_pos)
                                
                                # if it is not the first position of an exon and
                                # the current position is not the one right after/before the previous position,
                                # it means that the mutability for a given position(s) is missing
                                # then
                                # add a dictionary with all the alts and probability equals to 0,
                                # if there are more mutabilities of the consecutive positions missing, keep adding 0s

                                if pos != region[start]:
                                    expected_previous_pos = pos - 1
                                    # print(pos, prev_pos, expected_previous_pos)
                                    while prev_pos != expected_previous_pos:
                                        # print(pos, region[start], region[end], prev_pos, expected_previous_pos, cdna_pos)
                                        for altt in "ACGT":
                                            self.mutabilities_by_pos[cdna_pos][altt] = 0
                                        cdna_pos += update_pos
                                        expected_previous_pos -= 1

                                prev_pos = pos

                            # since at protein level we are looking at the nucleotide 
                            # changes of the translated codons we store them as they will be queried later
                            if self.reverse:
                                alt = transcribe[alt]

                            # add the mutability
                            self.mutabilities_by_pos[cdna_pos][alt] = mutability
                            # print(mutability, ref, alt, pos, cdna_pos, sep = "\t")

                        # segment_starting_pos = cdna_pos if not self.reverse else cdna_pos + segment_len
                        # print(self.chromosome, region[start], region[end], self.element, segment_len, cdna_pos, prev_pos, starting_cdna_pos, update_pos)

                        ##
                        #   IMPORTANT: filling step
                        ##
                        # check that all the positions at the end have been filled
                        # otherwise add 0s to them
                        if not self.reverse:
                            while cdna_pos < (starting_cdna_pos + segment_len):
                                for altt in "ACGT":
                                    self.mutabilities_by_pos[cdna_pos][altt] = 0
                                cdna_pos += 1
                        else:
                            while cdna_pos > (starting_cdna_pos - segment_len):
                                for altt in "ACGT":
                                    self.mutabilities_by_pos[cdna_pos][altt] = 0
                                cdna_pos -= 1

                        # print(self.chromosome, region[start], region[end], self.element, segment_len, cdna_pos, prev_pos, starting_cdna_pos, update_pos)
                        # print("\n")

                        # this is to get the cdna position pointer
                        # back to the biggest cdna position annotated so far
                        starting_cdna_pos = cdna_pos if not self.reverse else cdna_pos + segment_len


                    except ReaderError as e:
                        logger.warning(e.message)
                        continue
        except ReaderError as e:
            logger.warning("Reader error: %s. Regions being analysed %s", e.message, self.segments)
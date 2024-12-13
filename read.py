"""
Input file is a line-based format of records. First line of a record specifies
number of reads R and number haplotypes H to align against each read. Next R
lines specify reads, followed by H lines for haplotype sequences. The lines for
reads must contain five strings of equal length, separated by a single space.
The strings are: read sequence, read quality, insertion quality, deletion
quality, and gcp quality.

Example record:

2 2
CTGTGTCCATGTCAGAGCAATGGCCCAAGTCTGGGCCTGGG 888888888888'88888888'8888888888888'8'888 IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIN IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIN +++++++++++++++++++++++++++++++++++++++++
ATGTCAGAGCAATGGCCCAAGTCTGGGTCTGGG <AFFFKKKKKKKKKKKKKKKKKKKKKKKKKKKK IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIN IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIN +++++++++++++++++++++++++++++++++
CTGTGTCCATGTCAGAGCAACGGCCCAAGTCTGGGTCTGGG
CTGTGTCCATGTCAGAGCAATGGCCCAAGTCTGGGTCTGGG
"""

from typing import Union

class Read:
    def __init__(self, raw_data: Union[str, bytes]):
        num_columns = len(raw_data.split())
        if num_columns != 5:
            raise ValueError(f"Read data must contain 5 columns (read, read quality, insertion quality, deletion quality, gcp quality). Found {num_columns} columns.")

        # Check that all columns have the same length
        column_lengths = [len(column) for column in raw_data.split()]
        if len(set(column_lengths)) != 1:
            raise ValueError(f"All columns must have the same length.")

        self.raw_data = raw_data
        self.nucleotides = ''
        # These scores represent the probability of an error at each position in the read.
        # match to match ?
        self.base_qualities = []
        # These scores represent the probability of an insertion error at each position in the read.
        # match to insertion ?
        self.ins_qualities = []
        # These scores represent the probability of a deletion error at each position in the read.
        # match to deletion ?
        self.del_qualities = []
        # Gap-compressed probabilities ??
        # These qualities are used in bioinformatics to represent the
        # probability of gaps (insertions or deletions) in the sequence
        # alignment.
        # continue extending a gap ?
        self.gcp_qualities = []

        self.fill_nucleotides()
        self.fill_base_qualities()
        self.fill_ins_qualities()
        self.fill_del_qualities()
        self.fill_gcp_qualities()

    def fill_nucleotides(self):
        self.nucleotides = self.raw_data.split()[0]

    def fill_base_qualities(self):
        self.base_qualities = [self.phred_to_prob(ord(x)-33) for x in self.raw_data.split()[1]]

    def fill_ins_qualities(self):
        self.ins_qualities = [self.phred_to_prob(ord(x)-33) for x in self.raw_data.split()[2]]

    def fill_del_qualities(self):
        self.del_qualities = [self.phred_to_prob(ord(x)-33) for x in self.raw_data.split()[3]]

    def fill_gcp_qualities(self):
        self.gcp_qualities = [self.phred_to_prob(ord(x)-33) for x in self.raw_data.split()[4]]

    @staticmethod
    def phred_to_prob(phred: int) -> float:
        return 10 ** (-phred / 10)

    def __getitem__(self, i):
        return self.nucleotides[i]

    def __str__(self):
        s = f"Read({self.nucleotides.decode()})"
        s += f"\nBase qualities: {self.base_qualities}"
        s += f"\nInsertion qualities: {self.ins_qualities}"
        s += f"\nDeletion qualities: {self.del_qualities}"
        s += f"\nGCP qualities: {self.gcp_qualities}"
        return s

    def __repr__(self):
        return self.__str__()
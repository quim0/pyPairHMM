import matplotlib.pyplot as plt
import argparse

from filereader import PairHMMFile
from read import Read
import math

class PairHMMAligner:
    # std::numeric_limits<double>::max() / 16
    # INITIAL_CONSTANT = 1.7976931348623157e+308 / 16
    # std::numeric_limits<float>::max() / 16
    INITIAL_CONSTANT = 3.40282e+38 / 16
    def __init__(self, read: Read, haplotype: str):
        # Haplotide is in the x-axis, read is in the y-axis. (0,0) is the
        # top-left corner.
        self.read: Read = read
        self.haplotype: str = haplotype
        self.result = None

        self.M = [[0 for _ in range(len(haplotype)+1)] for _ in range(len(read.nucleotides)+1)]
        self.I = [[0 for _ in range(len(haplotype)+1)] for _ in range(len(read.nucleotides)+1)]
        self.D = [[0 for _ in range(len(haplotype)+1)] for _ in range(len(read.nucleotides)+1)]

    def prior(self, i, j):
        # Prior(i,j) = 1-Q_i if read[i] == haplotype[j] else Q_i/3
        if self.read[i-1] == self.haplotype[j-1]:
            return 1 - self.read.base_qualities[i-1]
        else:
            return self.read.base_qualities[i-1] / 3.0

    def align(self):
        if self.result is not None:
            return self.result

        # Initialize the first M row to 1/len(haplotype)
        # Other rows are initialized to 0
        for h in range(0, len(self.haplotype)):
            self.M[0][h] = 0.0
            self.I[0][h] = 0.0
            self.D[0][h] = self.INITIAL_CONSTANT / len(self.haplotype)

        # The first column is initialized to 0 for all matrices
        # Start at row 1 because row 0 is already initialized
        for r in range(1, len(self.read.nucleotides)):
            self.M[r][0] = 0.0
            self.I[r][0] = 0.0
            self.D[r][0] = 0.0


        for r in range(1, len(self.read.nucleotides)+1):
            for h in range(1, len(self.haplotype)+1):
                # Define the transition probabilities
                match_continuation_prob = 1.0 - (self.read.ins_qualities[r-1] + self.read.del_qualities[r-1]) # alpha
                gap_to_match_prob = 1.0 - self.read.gcp_qualities[r-1] # beta
                match_to_insertion_prob = self.read.ins_qualities[r-1]
                # Why the if? It's written like that in the pairhmm code
                match_to_deletion_prob = self.read.del_qualities[r-1] if r < len(self.read.nucleotides) else 1.0
                gap_continuation_prob = self.read.gcp_qualities[r-1] if r < len(self.read.nucleotides) else 1.0

                print(f'Probabilities:')
                print(f'match_continuation_prob: {match_continuation_prob}')
                print(f'gap_to_match_prob: {gap_to_match_prob}')
                print(f'match_to_insertion_prob: {match_to_insertion_prob}')
                print(f'match_to_deletion_prob: {match_to_deletion_prob}')
                print(f'gap_continuation_prob: {gap_continuation_prob}')
                print(f'Prior: {self.prior(r, h)} match? {self.read[r-1]} == {self.haplotype[h-1]}')

                self.M[r][h] = self.prior(r, h) * (
                    match_continuation_prob * self.M[r-1][h-1] +
                    gap_to_match_prob * (self.I[r-1][h-1] + self.D[r-1][h-1])
                )
                self.I[r][h] = match_to_insertion_prob * self.M[r-1][h] + gap_continuation_prob * self.I[r-1][h]
                self.D[r][h] = match_to_deletion_prob * self.M[r][h-1] + gap_continuation_prob * self.D[r][h-1]

                #if self.M[r][h] < 1e29:
                #    self.M[r][h] = 0
                #if self.I[r][h] < 1e29:
                #    self.I[r][h] = 0
                #if self.D[r][h] < 1e29:
                #    self.D[r][h] = 0

        # The result is the sum of the last row from the M and I matrices
        res = 0
        for h in range(0, len(self.haplotype)+1):
            res += self.M[len(self.read.nucleotides)][h] + self.I[len(self.read.nucleotides)][h]

        res = math.log10(res) - math.log10(self.INITIAL_CONSTANT)

        self.result = res
        return self.result

    def plot(self):
        # One axes for each matrix
        fig, ax = plt.subplots(3, 1, figsize=(15, 22))

        ax[0].imshow(self.M, cmap='hot')
        ax[0].set_title('M matrix')

        ax[1].imshow(self.I, cmap='hot')
        ax[1].set_title('I matrix')

        ax[2].imshow(self.D, cmap='hot')
        ax[2].set_title('D matrix')

        plt.show()

if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description="PairHMM")
    arg_parser.add_argument("file_path", type=str, help="Path to the input file")
    args = arg_parser.parse_args()

    pair_hmm_file = PairHMMFile(args.file_path)
    for batch in pair_hmm_file.batchs:
        for read in batch.reads:
            for haplotype in batch.haplotypes:
                print(f'read: {read.nucleotides}')
                aligner = PairHMMAligner(read, haplotype)
                print(f'read length: {len(aligner.read.nucleotides)}, haplotype length: {len(aligner.haplotype)}')
                print(aligner.align())
                aligner.plot()
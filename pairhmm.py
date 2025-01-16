import matplotlib.pyplot as plt
import numpy as np
import argparse
import re

from filereader import PairHMMFile, get_file_data
from z3 import *
from read import Read
import math

def z3_output_to_latex(expression):
    """
    Convert a Z3-like mathematical expression to LaTeX format,
    including Greek symbols and subscripts.
    """
    # Replace subscripts with LaTeX-style subscripts
    expression = re.sub(r'(\w+)_(\d+,\d+)', r'\1_{\2}', expression)
    expression = re.sub(r'(\w+)_(\d+)', r'\1_{\2}', expression)

    # Replace multiplication (*) with LaTeX multiplication (\cdot)
    expression = expression.replace('*', r' \cdot ')

    # Replace addition (+) with spaces around it for better formatting
    expression = expression.replace('+', r' + ')

    # Replace specific terms with LaTeX Greek symbols
    expression = expression.replace('alpha', r'\alpha')
    expression = expression.replace('beta', r'\beta')
    expression = expression.replace('delta', r'\delta')
    expression = expression.replace('epsilon', r'\epsilon')
    expression = expression.replace('prior', r'\lambda')  # Map prior to \lambda

    # Wrap the entire expression in math mode
    return f"\\[ {expression} \\]"

class PairHMMAligner:
    # std::numeric_limits<float>::max() / 16
    INITIAL_CONSTANT_FP32 = 3.40282e+38 / 16
    # std::numeric_limits<double>::max() / 16
    INITIAL_CONSTANT_FP64 = 1.7976931348623157e+308 / 16
    def __init__(self, read: Read, haplotype: str, precision: str = 'float'):
        # Haplotide is in the x-axis, read is in the y-axis. (0,0) is the
        # top-left corner.
        self.read: Read = read
        self.haplotype: str = haplotype
        self.result = None

        if precision not in ['float', 'double']:
            raise ValueError('Precision must be float or double')

        self.INITIAL_CONSTANT = self.INITIAL_CONSTANT_FP32 if precision == 'float' else self.INITIAL_CONSTANT_FP64

        arr_dtype = 'd' if precision == 'double' else 'f'
        self.M = np.array([[0.0 for _ in range(len(haplotype)+1)] for _ in range(len(read.nucleotides)+1)], dtype=arr_dtype)
        self.I = np.array([[0.0 for _ in range(len(haplotype)+1)] for _ in range(len(read.nucleotides)+1)], dtype=arr_dtype)
        self.D = np.array([[0.0 for _ in range(len(haplotype)+1)] for _ in range(len(read.nucleotides)+1)], dtype=arr_dtype)

    def prior(self, i, j):
        # Prior(i,j) = 1-Q_i if read[i] == haplotype[j] else Q_i/3
        if self.read[i-1] == self.haplotype[j-1]:
            return 1 - self.read.base_qualities[i-1]
        else:
            return self.read.base_qualities[i-1] / 3.0

    def initialize_matrices(self):
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

    def align(self):
        if self.result is not None:
            return self.result

        self.initialize_matrices()

        for r in range(1, len(self.read.nucleotides)+1):
            for h in range(1, len(self.haplotype)+1):
                # Define the transition probabilities
                match_continuation_prob = 1.0 - (self.read.ins_qualities[r-1] + self.read.del_qualities[r-1]) # alpha
                gap_to_match_prob = 1.0 - self.read.gcp_qualities[r-1] # beta
                match_to_insertion_prob = self.read.ins_qualities[r-1]
                # Why the if? It's written like that in the pairhmm code
                match_to_deletion_prob = self.read.del_qualities[r-1] if r < len(self.read.nucleotides) else 1.0
                gap_continuation_prob = self.read.gcp_qualities[r-1] if r < len(self.read.nucleotides) else 1.0

                #print(f'Probabilities:')
                #print(f'match_continuation_prob: {match_continuation_prob}')
                #print(f'gap_to_match_prob: {gap_to_match_prob}')
                #print(f'match_to_insertion_prob: {match_to_insertion_prob}')
                #print(f'match_to_deletion_prob: {match_to_deletion_prob}')
                #print(f'gap_continuation_prob: {gap_continuation_prob}')
                #print(f'Prior: {self.prior(r, h)} match? {self.read[r-1]} == {self.haplotype[h-1]}')

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
        fig, ax = plt.subplots(1, 3, figsize=(20, 30))

        #ax[0].imshow(np.log10(self.M) , cmap='hot')
        ax[0].imshow(self.M, cmap='hot')
        ax[0].set_title('M matrix')

        # highlight the maximum value in each row
        for r in range(0, len(self.read.nucleotides)+1):
            max_idx = np.argmax(self.M[r])
            ax[0].plot(max_idx, r, 'ro')

        # Put the letters in the ticks
        # 1. Set the origin on the top-left corner
        # set x ticks in the top axis
        ax[0].xaxis.tick_top()
        # 2. Set the ticks
        ax[0].set_xticks(np.arange(len(self.haplotype)))
        ax[0].set_xticklabels(list(self.haplotype))
        ax[0].set_yticks(np.arange(len(self.read.nucleotides)))
        ax[0].set_yticklabels(list(self.read.nucleotides))

        ax[1].imshow(np.log10(self.I) , cmap='hot')
        #ax[1].imshow(self.I, cmap='hot')
        ax[1].set_title('I matrix')

        ax[2].imshow(np.log10(self.D) , cmap='hot')
        #ax[2].imshow(self.D, cmap='hot')
        ax[2].set_title('D matrix')

        plt.show()

    def _is_monotone(self, matrix):
        # For each row, get the index of the maximum value, to be monotone the matrix must:
        # - The index of a row must be less or equal than the index of the next row
        for ridx, r in enumerate(matrix):
            max_idx = np.argmax(r)
            if ridx > 0 and max_idx < np.argmax(matrix[ridx-1]):
                return False
        return True

    def is_monotone(self):
        return (self._is_monotone(self.M), self._is_monotone(self.I), self._is_monotone(self.D))

class SymbolicPairHMMAligner(PairHMMAligner):
    def __init__(self, read: Read, haplotype: str):
        super().__init__(read, haplotype)

    def prior(self, i, j):
        # Prior(i,j) = 1-Q_i if read[i] == haplotype[j] else Q_i/3
        if self.read[i-1] == self.haplotype[j-1]:
            return 1 - self.read.base_qualities[i-1]
        else:
            return self.read.base_qualities[i-1] / 3.0

    def initialize_matrices(self):
        for h in range(0, len(self.haplotype)):
            self.M[0][h] = 0.0
            self.I[0][h] = 0.0
            self.D[0][h] = Real(f'initial')


    def align(self):
        if self.result is not None:
            return self.result

        # Create the variables
        self.M = [[0.0 for j in range(len(self.haplotype)+1)] for i in range(len(self.read.nucleotides)+1)]
        self.I = [[0.0 for j in range(len(self.haplotype)+1)] for i in range(len(self.read.nucleotides)+1)]
        self.D = [[0.0 for j in range(len(self.haplotype)+1)] for i in range(len(self.read.nucleotides)+1)]

        self.initialize_matrices()

        # Add the rest of the constraints
        for r in range(1, len(self.read.nucleotides)+1):
            for h in range(1, len(self.haplotype)+1):
                # Define the transition probabilities
                # match_continuation_prob = 1.0 - (2.0 * Real(f'delta_{r}')) # alpha
                match_continuation_prob = Real(f'alpha_{r}')
                #gap_to_match_prob = 1.0 - Real('epsilon') # beta
                gap_to_match_prob = Real(f'beta')
                match_to_insertion_prob = Real('epsilon')
                # Why the if? It's written like that in the pairhmm code
                match_to_deletion_prob = Real(f'delta_{r}') if r < len(self.read.nucleotides) else 1.0

                gap_continuation_prob = Real('epsilon') if r < len(self.read.nucleotides) else 1.0

                prior = Real(f'prior_{r},{h}')

                self.M[r][h] = prior * (
                    match_continuation_prob * self.M[r-1][h-1] +
                    gap_to_match_prob * (self.I[r-1][h-1] + self.D[r-1][h-1])
                )

                self.I[r][h] = match_to_insertion_prob * self.M[r-1][h] + gap_continuation_prob * self.I[r-1][h]

                self.D[r][h] = match_to_deletion_prob * self.M[r][h-1] + gap_continuation_prob * self.D[r][h-1]


# Create subclass of pairHmmaligner
class PairHMMHeuristicAligner(PairHMMAligner):
    def __init__(self, read: Read, haplotype: str, precision: str = 'float'):
        super().__init__(read, haplotype, precision)

    def align(self):
        self.initialize_matrices()

        # Iterate over all the diagonals
        # Number of diagonals = haplotype length
        for d2 in range(0, len(self.haplotype)+2):
            d = len(self.haplotype)+1 - d2
            v = 0
            h = d
            h += 1
            v += 1
            while h < len(self.haplotype)+1 and v < len(self.read.nucleotides)+1:
                #and self.haplotype[h-1] == self.read[v-1]:

                r = v
                h = h

                # Define the transition probabilities
                match_continuation_prob = 1.0 - (self.read.ins_qualities[r-1] + self.read.del_qualities[r-1]) # alpha
                gap_to_match_prob = 1.0 - self.read.gcp_qualities[r-1] # beta
                match_to_insertion_prob = self.read.ins_qualities[r-1]
                # Why the if? It's written like that in the pairhmm code
                match_to_deletion_prob = self.read.del_qualities[r-1] if r < len(self.read.nucleotides) else 1.0
                gap_continuation_prob = self.read.gcp_qualities[r-1] if r < len(self.read.nucleotides) else 1.0

                self.M[r][h] = self.prior(r, h) * (
                    match_continuation_prob * self.M[r-1][h-1] +
                    gap_to_match_prob * (self.I[r-1][h-1] + self.D[r-1][h-1])
                )
                self.I[r][h] = match_to_insertion_prob * self.M[r-1][h] + gap_continuation_prob * self.I[r-1][h]
                self.D[r][h] = match_to_deletion_prob * self.M[r][h-1] + gap_continuation_prob * self.D[r][h-1]

                if self.M[r][h] < 1e29:
                    break

                v += 1
                h += 1

        # Replace all the zeros in the last row of the M and I matrix with 1e20
        self.M[self.M == 0]= 1e20
        self.I[self.I == 0]= 1e20
        self.D[self.D == 0]= 1e20

        # The result is the sum of the last row from the M and I matrices
        res = 0
        for h in range(0, len(self.haplotype)+1):
            res += self.M[len(self.read.nucleotides)][h] + self.I[len(self.read.nucleotides)][h]

        res = math.log10(res) - math.log10(self.INITIAL_CONSTANT)

        self.result = res
        return self.result

if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(description="PairHMM")
    arg_parser.add_argument("file_path", type=str, help="Path to the input file")
    arg_parser.add_argument("--precision", type=str, default="float", help="Precision to use (float or double)")
    args = arg_parser.parse_args()


    #print(
    #    get_file_data(args.file_path)
    #)

    #exit()

    pair_hmm_file = PairHMMFile(args.file_path, num_batchs=10)
    print(f'Max read length: {pair_hmm_file.max_read_length()}')
    print(f'Max haplotype length: {pair_hmm_file.max_haplotype_length()}')

    #exit()

    for batch in pair_hmm_file.batchs:
        for read in batch.reads:
            for haplotype in batch.haplotypes:
                print(f'read: {read.nucleotides}')
                aligner = PairHMMAligner(read, haplotype)
                #aligner = SymbolicPairHMMAligner(read, haplotype)
                #heuristic_aligner = PairHMMHeuristicAligner(read, haplotype)
                exact_res = aligner.align()
                #heuristic_res = heuristic_aligner.align()
                print(f'read length: {len(aligner.read.nucleotides)}, haplotype length: {len(aligner.haplotype)}')
                print(f'Exact: {aligner.align()}')

                # get result only considering the maximum  value of the last row
                print(f'Exact2: {math.log10(np.max(aligner.M[-1])) - math.log10(aligner.INITIAL_CONSTANT)}')

                #print(z3_output_to_latex(str(simplify(aligner.M[5][5]))))
                #exit(0)
                #print(f'Heuristic: {heuristic_aligner.align()}')
                #print(f'Difference: {exact_res - heuristic_res}')
                #aligner.plot()
                #heuristic_aligner.plot()
                #print(f'Is monotone: {aligner.is_monotone()}')
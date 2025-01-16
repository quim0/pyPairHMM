from read import Read

class PairHMMBatch:
    def __init__(self, reads: list[Read], haplotypes: list[str]):
        self.reads: list[Read] = reads
        self.haplotypes: list[str] = haplotypes

    def max_read_length(self):
        return max([len(read.nucleotides) for read in self.reads])

    def max_haplotype_length(self):
        return max([len(haplotype) for haplotype in self.haplotypes])

    def __str__(self):
        return f"PairHMMBatch(num_reads={len(self.reads)}, num_haplotypes={len(self.haplotypes)})"

    def __repr__(self):
        return str(self)

class PairHMMFile:
    def __init__(self, file_path: str, num_batchs: int = 0):
        self.batchs: list[PairHMMBatch] = []
        self.file_path: str = file_path
        self.num_batchs_to_read: int = num_batchs
        self.current = 0

        self.parse_file()

    def max_read_length(self):
        return max([batch.max_read_length() for batch in self.batchs])

    def max_haplotype_length(self):
        return max([batch.max_haplotype_length() for batch in self.batchs])

    def parse_file(self):
        with open(self.file_path, 'r') as f:
            lines = f.readlines()

            i = 0

            while i < len(lines):
                curr_reads = []
                curr_haplotypes = []

                # go to the first non-empty line
                while i < len(lines) and lines[i].strip() == '':
                    i += 1

                # read the number of reads and haplotypes
                num_reads, num_haplotypes = map(int, lines[i].split())

                # read the reads
                i += 1
                for _ in range(num_reads):
                    read = Read(lines[i])
                    curr_reads.append(read)
                    i += 1

                # read the haplotypes
                for _ in range(num_haplotypes):
                    haplotype = lines[i].strip()
                    curr_haplotypes.append(haplotype)
                    i += 1

                self.batchs.append(PairHMMBatch(curr_reads, curr_haplotypes))

                if len(self.batchs) == self.num_batchs_to_read:
                    break

    def __str__(self):
        return f"PairHMMFile({self.file_path}, num_batchs={len(self.batchs)})"

    def __repr__(self):
        return str(self)

def get_file_data(file_path: str) -> dict:
    num_batchs = 0
    max_read_length = 0
    max_haplotype_length = 0

    max_num_reads = 0
    max_num_haplotypes = 0

    with open(file_path, 'r') as f:
        for line in f:
            if line.strip() == '':
                continue

            num_batchs += 1

            num_reads, num_haplotypes = map(int, line.split())

            max_num_reads = max(max_num_reads, num_reads)
            max_num_haplotypes = max(max_num_haplotypes, num_haplotypes)

            for _ in range(num_reads):
                line = next(f)
                read = line.strip()

                max_read_length = max(max_read_length, len(read))

            for _ in range(num_haplotypes):
                line = next(f)
                haplotype = line.strip()

                max_haplotype_length = max(max_haplotype_length, len(haplotype))

    return {
        'num_batchs': num_batchs,
        'max_read_length': max_read_length,
        'max_haplotype_length': max_haplotype_length,
        'max_num_reads': max_num_reads,
        'max_num_haplotypes': max_num_haplotypes
    }
from read import Read

class PairHMMBatch:
    def __init__(self, reads: list[Read], haplotypes: list[str]):
        self.reads: list[Read] = reads
        self.haplotypes: list[str] = haplotypes

    def __str__(self):
        return f"PairHMMBatch(num_reads={len(self.reads)}, num_haplotypes={len(self.haplotypes)})"

    def __repr__(self):
        return str(self)

class PairHMMFile:
    def __init__(self, file_path: str):
        self.batchs: list[PairHMMBatch] = []
        self.file_path: str = file_path
        self.current = 0

        self.parse_file()

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

    def __str__(self):
        return f"PairHMMFile({self.file_path}, num_batchs={len(self.batchs)})"

    def __repr__(self):
        return str(self)
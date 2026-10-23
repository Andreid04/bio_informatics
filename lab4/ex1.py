import random
from typing import List, Tuple
#./GCF_001399695.1_ASM139969v1_cds_from_genomic.fna

def load_fna_genome(filepath: str) -> str:
    chunks = []
    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            chunks.append(line.upper())
    return "".join(chunks)


def sample_genome_uniform(
    genome: str,
    num_samples: int = 20000,
    min_len: int = 100,
    max_len: int = 150,
) -> List[Tuple[str, int, int]]:
    total_len = len(genome)
    max_valid_start = total_len - min_len

    if max_valid_start <= 0:
        raise ValueError(f"Genome length shorter than min_len.")

    samples = []
    for _ in range(num_samples):
        frag_len = random.randint(min_len, max_len)
        valid_end_start = total_len - frag_len

        if valid_end_start < 0:
            frag_len = total_len
            start_idx = 0
        else:
            start_idx = random.randint(0, valid_end_start)

        seq = genome[start_idx : start_idx + frag_len]
        samples.append((seq, start_idx, frag_len))

    return samples

def assemble_simple(reads: List[str], min_overlap: int) -> str:
    lookup = {}
    for r in reads:
        lookup[r[:min_overlap]] = r

    current = reads[0]
    assembled = current
    used = {current}

    while True:
        suffix = current[-min_overlap:]
        nxt = lookup.get(suffix)
        if nxt and nxt not in used:
            assembled += nxt[min_overlap:]
            used.add(nxt)
            current = nxt
        else:
            break

    return assembled


def main():
    input_file = "./GCF_001399695.1_ASM139969v1_cds_from_genomic.fna"
    num_samples = 20000
    min_len = 100
    max_len = 150

    genome = load_fna_genome(input_file)
    reads_vector: List[Tuple[str, int, int]] = sample_genome_uniform(
        genome,
        num_samples=num_samples,
        min_len=min_len,
        max_len=max_len,
    )

    blind_reads: List[str] = [seq for seq, _, _ in reads_vector]
    random.shuffle(blind_reads)

    reconstructed = assemble_simple(blind_reads, min_overlap = 8)

    print(f"Total genome length: {len(genome)}")
    print(f"Total samples in vector: {len(reads_vector)}")
    print(f"Length range: [{min_len}, {max_len}]")
    print(f"Sample 1 length: {reads_vector[0][2]}, start: {reads_vector[0][1]}")
    print(f"Sample {len(reads_vector)} length: {reads_vector[-1][2]}, start: {reads_vector[-1][1]}")
    print(f"Reconstructed sequence length: {len(reconstructed)}")
    print(f"the reconstructed genome is: \n {reconstructed}")

if __name__ == "__main__":
    main()
import math
import random
from typing import List, Tuple
import urllib.request


def fetch_ncbi_sequence(accession: str) -> str:
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req) as resp:
        lines = resp.read().decode("utf-8").splitlines()
    return "".join(line.strip().upper() for line in lines if not line.startswith(">"))


def sample_fragments(
    dna: str,
    num_samples: int = 10,
    min_len: int = 100,
    max_len: int = 3000,
) -> List[Tuple[int, str]]:
    total_len = len(dna)
    samples = []
    for _ in range(num_samples):
        target_len = random.randint(min_len, min(max_len, total_len))
        start_idx = random.randint(0, total_len - target_len)
        fragment = dna[start_idx : start_idx + target_len]
        samples.append((len(fragment), fragment))
    return samples


def simulate_gel_migration(
    fragment_lengths: List[int],
    gel_length_cm: float = 10.0,
) -> List[Tuple[int, float]]:
    min_bp = 100
    max_bp = 3000
    log_min = math.log10(min_bp)
    log_max = math.log10(max_bp)

    results = []
    for bp in fragment_lengths:
        log_bp = math.log10(bp)
        normalized = (log_max - log_bp) / (log_max - log_min)
        distance = round(1.0 + normalized * (gel_length_cm - 2.0), 2)
        results.append((bp, distance))

    results.sort(key=lambda x: x[1])
    return results


def main():
    accession = "NM_000518.5"
    try:
        dna = fetch_ncbi_sequence(accession)
    except Exception:
        random.seed(42)
        dna = "".join(random.choices(["A", "C", "G", "T"], k=2000))

    if len(dna) > 3000:
        dna = dna[:3000]

    samples = sample_fragments(dna, num_samples=10, min_len=100, max_len=3000)
    lengths_array = [bp for bp, _ in samples]

    simulation_results = simulate_gel_migration(lengths_array)

    print(f"Original Sequence Length: {len(dna)} bp")
    print(f"Sampled Fragments Count: {len(samples)}")
    print(f"Fragment Sizes (bp): {lengths_array}")
    print("Gel Electrophoresis Simulation (Top -> Bottom):")
    print(f"{'Well / Cathode (-)':^35}")
    print("-" * 35)
    for bp, dist in simulation_results:
        marker = " " * int(dist * 3) + f"■ ({bp} bp, {dist} cm)"
        print(marker)
    print("-" * 35)
    print(f"{'Anode (+)':^35}")


if __name__ == "__main__":
    main()
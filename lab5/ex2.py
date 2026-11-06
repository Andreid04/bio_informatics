import random
import re
from typing import Dict, List, Tuple
import urllib.request


def fetch_ncbi_sequence(accession: str) -> str:
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req) as resp:
        lines = resp.read().decode("utf-8").splitlines()
    return "".join(line.strip().upper() for line in lines if not line.startswith(">"))


def detect_tandem_repeats(
    dna: str,
    min_unit: int = 3,
    max_unit: int = 5,
    min_total: int = 6,
    max_total: int = 10,
) -> List[Dict]:
    results = []
    n = len(dna)

    for unit_len in range(min_unit, max_unit + 1):
        for start in range(n - unit_len):
            unit = dna[start : start + unit_len]
            k = 1
            while start + (k + 1) * unit_len <= n and dna[start + k * unit_len : start + (k + 1) * unit_len] == unit:
                k += 1

            total_span = k * unit_len
            if min_total <= total_span <= max_total and k >= 2:
                results.append(
                    {
                        "unit": unit,
                        "unit_length": unit_len,
                        "repeat_count": k,
                        "total_length": total_span,
                        "start_pos": start,
                        "end_pos": start + total_span,
                        "full_match": unit * k,
                    }
                )

    dedup = {}
    for r in results:
        key = (r["start_pos"], r["end_pos"], r["unit"])
        dedup[key] = r

    unique_results = sorted(dedup.values(), key=lambda x: x["start_pos"])
    return unique_results


def main():
    accession = "NM_000518.5"
    try:
        dna = fetch_ncbi_sequence(accession)
    except Exception:
        random.seed(123)
        dna = "".join(random.choices(["A", "C", "G", "T"], k=2000))

    if len(dna) > 3000:
        dna = dna[:3000]

    repeats = detect_tandem_repeats(
        dna, min_unit=3, max_unit=5, min_total=6, max_total=10
    )

    print(f"Sequence Length: {len(dna)} bp")
    print(f"Total Repeated Units Detected: {len(repeats)}")
    print(f"{'Start':<8} {'End':<8} {'Unit':<8} {'UnitLen':<10} {'Copies':<8} {'TotalLen':<10} {'Sequence'}")
    print("-" * 65)
    for r in repeats:
        print(
            f"{r['start_pos']:<8} {r['end_pos']:<8} {r['unit']:<8} {r['unit_length']:<10} {r['repeat_count']:<8} {r['total_length']:<10} {r['full_match']}"
        )


if __name__ == "__main__":
    main()
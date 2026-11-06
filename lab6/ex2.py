from collections import defaultdict
from bisect import bisect_left

def read_fasta(filename):
    sequence = []
    header = ""

    with open(filename, "r") as file:
        for line in file:
            line = line.strip()

            if line.startswith(">"):
                if not header:
                    header = line
                elif sequence:
                    raise ValueError(
                        "Process each chromosome separately."
                    )
            else:
                sequence.append(line.upper())

    genome = "".join(sequence)

    return header, genome


def reverse_complement(sequence):
    complement = str.maketrans("ACGT", "TGCA")
    return sequence.translate(complement)[::-1]


def detect_inverted_repeats(
    genome,
    window_size=15,
    min_distance=100
):

    detections = []

    genome_length = len(genome)

    positions = defaultdict(list)

    for j in range(genome_length - window_size + 1):
        s2 = genome[j:j + window_size]

        if all(base in "ACGT" for base in s2):
            positions[s2].append(j)

    for i in range(genome_length - window_size + 1):

        s1 = genome[i:i + window_size]

        if not all(base in "ACGT" for base in s1):
            continue

        target = reverse_complement(s1)

        minimum_s2_position = i + window_size + min_distance

        matching_positions = positions.get(target, [])

        index = bisect_left(
            matching_positions,
            minimum_s2_position
        )

        for j in matching_positions[index:]:

            detections.append({
                "S1_start": i,
                "S1_end": i + window_size,
                "S1_sequence": s1,
                "S2_start": j,
                "S2_end": j + window_size,
                "S2_sequence": target,
                "distance": j - (i + window_size)
            })

    return detections


def main():

    print("INVERTED REPEAT DETECTION")
    print("=" * 60)

    filename = "./GCF_001399695.1_ASM139969v1_cds_from_genomic.fna"

    header, genome = read_fasta(filename)

    print("\nGenome:", header)
    print("Genome length:", len(genome), "bases")

    window_size = 15
    min_distance = 100

    print("\nParameters:")
    print("S1/S2 window size:", window_size)
    print("Minimum distance:", min_distance)

    print("\nSearching genome...")

    results = detect_inverted_repeats(
        genome,
        window_size,
        min_distance
    )

    print("\nDetection completed.")
    print("Total candidate pairs:", len(results))

    print("\nFirst 20 results:")
    print("-" * 100)

    print(
        f"{'S1 Start':<12}"
        f"{'S1 End':<12}"
        f"{'S2 Start':<12}"
        f"{'S2 End':<12}"
        f"{'Distance':<12}"
    )

    for result in results[:20]:

        print(
            f"{result['S1_start']:<12}"
            f"{result['S1_end']:<12}"
            f"{result['S2_start']:<12}"
            f"{result['S2_end']:<12}"
            f"{result['distance']:<12}"
        )

    # Save all detections to a CSV file.
    with open("inverted_repeats.csv", "w") as file:

        file.write(
            "S1_start,S1_end,S1_sequence,"
            "S2_start,S2_end,S2_sequence,distance\n"
        )

        for result in results:

            file.write(
                f"{result['S1_start']},"
                f"{result['S1_end']},"
                f"{result['S1_sequence']},"
                f"{result['S2_start']},"
                f"{result['S2_end']},"
                f"{result['S2_sequence']},"
                f"{result['distance']}\n"
            )

    print("\nResults saved to inverted_repeats.csv")


if __name__ == "__main__":
    main()
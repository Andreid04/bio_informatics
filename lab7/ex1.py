
import csv

ENZYMES = {
    "EcoRI": {
        "sequence": "GAATTC",
        "cut": 1
    },
    "BamHI": {
        "sequence": "GGATCC",
        "cut": 1
    },
    "HindIII": {
        "sequence": "AAGCTT",
        "cut": 1
    },
    "TagI": {
        "sequence": "TCGA",
        "cut": 1
    },
    "HaeIII": {
        "sequence": "GGCC",
        "cut": 2
    }
}


FASTA_FILE = "./GCF_001399695.1_ASM139969v1_cds_from_genomic.fna"


def read_fasta(filename):

    records = []
    header = None
    sequence = []

    with open(filename, "r") as file:

        for line in file:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):

                if header is not None:
                    records.append(
                        (header, "".join(sequence))
                    )

                header = line[1:]
                sequence = []

            else:
                sequence.append(line.upper())

    if header is not None:
        records.append(
            (header, "".join(sequence))
        )

    return records


def find_cleavage_sites(sequence, selected_enzymes):

    cuts = {}

    enzyme_hits = []

    for enzyme_name in selected_enzymes:

        enzyme = ENZYMES[enzyme_name]

        recognition = enzyme["sequence"]
        offset = enzyme["cut"]

        start = 0

        while True:

            position = sequence.find(recognition, start)

            if position == -1:
                break

            cleavage_position = position + offset

            cuts.setdefault(cleavage_position, []).append(
                enzyme_name
            )

            enzyme_hits.append({
                "enzyme": enzyme_name,
                "recognition": recognition,
                "site_start": position,
                "cut_position": cleavage_position
            })

            # allow overlapping matches
            start = position + 1

    return cuts, enzyme_hits


def calculate_fragments(sequence_length, cuts):

    cut_positions = sorted(cuts.keys())

    boundaries = [0] + cut_positions + [sequence_length]

    fragments = []

    for i in range(len(boundaries) - 1):

        start = boundaries[i]
        end = boundaries[i + 1]

        fragments.append({
            "fragment": i + 1,
            "start": start,
            "end": end,
            "length": end - start
        })

    return fragments


def digest_sequence(sequence, selected_enzymes):

    cuts, enzyme_hits = find_cleavage_sites(
        sequence,
        selected_enzymes
    )

    fragments = calculate_fragments(
        len(sequence),
        cuts
    )

    return cuts, enzyme_hits, fragments


def save_results(all_cleavages, all_fragments):

    with open("restriction_cleavages.csv", "w", newline="") as file:

        writer = csv.writer(file)

        writer.writerow([
            "Sequence",
            "Enzyme",
            "Recognition",
            "Site_start_1based",
            "Cut_after_base"
        ])

        for row in all_cleavages:
            writer.writerow(row)

    with open("restriction_fragments.csv", "w", newline="") as file:

        writer = csv.writer(file)

        writer.writerow([
            "Sequence",
            "Fragment",
            "Start_1based",
            "End_1based",
            "Length"
        ])

        for row in all_fragments:
            writer.writerow(row)


def main():

    print("-" * 10)
    print("rest enzyme digestion")
    print("-" * 10)

    records = read_fasta(FASTA_FILE)

    print("Number of sequences:", len(records))

    print("\nAvailable restriction enzymes:")

    for name, enzyme in ENZYMES.items():

        print(
            f"{name:<10} "
            f"Recognition: {enzyme['sequence']:<8} "
            f"Cut offset: {enzyme['cut']}"
        )

    print("\nEnter enzyme names separated by commas.")
    print("Example: EcoRI,BamHI,HindIII")

    user_input = input("\nSelected enzymes: ")

    selected_enzymes = [
        name.strip()
        for name in user_input.split(",")
        if name.strip()
    ]

    selected_enzymes = list(dict.fromkeys(selected_enzymes))

    invalid = [
        name for name in selected_enzymes
        if name not in ENZYMES
    ]

    if invalid:
        print("Invalid enzyme names:", invalid)
        return

    if not selected_enzymes:
        print("No enzymes selected.")
        return

    print("\nSelected enzymes:", ", ".join(selected_enzymes))

    all_cleavages = []
    all_fragments = []

    total_cleavage_events = 0
    total_unique_cuts = 0

    for record_name, sequence in records:

        cuts, hits, fragments = digest_sequence(
            sequence,
            selected_enzymes
        )

        total_cleavage_events += len(hits)
        total_unique_cuts += len(cuts)

        print("\n" + "-" * 10)
        print("Sequence:", record_name[:70])
        print("Sequence length:", len(sequence), "bp")

        print("\nCLEAVAGE RESULTS")
        print("-" * 10)

        print("Total enzyme recognition events:", len(hits))
        print("Unique cleavage positions:", len(cuts))

        if not hits:
            print("No restriction sites found.")

        else:

            print(
                f"\n{'Enzyme':<12}"
                f"{'Site start':<15}"
                f"{'Cut after':<15}"
                f"{'Recognition':<15}"
            )

            for hit in hits:

                print(
                    f"{hit['enzyme']:<12}"
                    f"{hit['site_start'] + 1:<15}"
                    f"{hit['cut_position']:<15}"
                    f"{hit['recognition']:<15}"
                )

                all_cleavages.append([
                    record_name,
                    hit["enzyme"],
                    hit["recognition"],
                    hit["site_start"] + 1,
                    hit["cut_position"]
                ])

            print("\nCOMBINED CUT POSITIONS")
            print("-" * 10)

            for position, enzymes in sorted(cuts.items()):

                print(
                    f"Position {position}: "
                    f"{', '.join(enzymes)}"
                )

        print("\nDNA FRAGMENTS")
        print("-" * 10)

        print(
            f"{'Fragment':<12}"
            f"{'Start':<12}"
            f"{'End':<12}"
            f"{'Length (bp)':<15}"
        )

        for fragment in fragments:

            print(
                f"{fragment['fragment']:<12}"
                f"{fragment['start'] + 1:<12}"
                f"{fragment['end']:<12}"
                f"{fragment['length']:<15}"
            )

            all_fragments.append([
                record_name,
                fragment["fragment"],
                fragment["start"] + 1,
                fragment["end"],
                fragment["length"]
            ])

        # Verify that fragment lengths add up to the DNA length
        assert sum(f["length"] for f in fragments) == len(sequence)

    print("\n" + "-" * 10)
    print("Final Summary")
    print("-" * 10)

    print("Sequences analyzed:", len(records))
    print("Total recognition events:", total_cleavage_events)
    print("Total unique cleavage positions:", total_unique_cuts)

    save_results(all_cleavages, all_fragments)

    print("\nResults exported:")
    print("- restriction_cleavages.csv")
    print("- restriction_fragments.csv")


if __name__ == "__main__":
    main()

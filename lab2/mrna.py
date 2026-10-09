

mrna = "AAUGGCUUUUUGUUGUGUUGCUUAG"

genetic_code = {
    # UU
    "UUU": "Phe", "UUC": "Phe",
    "UUA": "Leu", "UUG": "Leu",

    # UC
    "UCU": "Ser", "UCC": "Ser",
    "UCA": "Ser", "UCG": "Ser",

    # UA
    "UAU": "Tyr", "UAC": "Tyr",
    "UAA": "Stop", "UAG": "Stop",

    # UG
    "UGU": "Cys", "UGC": "Cys",
    "UGA": "Stop", "UGG": "Trp",

    # CU
    "CUU": "Leu", "CUC": "Leu",
    "CUA": "Leu", "CUG": "Leu",

    # CC
    "CCU": "Pro", "CCC": "Pro",
    "CCA": "Pro", "CCG": "Pro",

    # CA
    "CAU": "His", "CAC": "His",
    "CAA": "Gln", "CAG": "Gln",

    # CG
    "CGU": "Arg", "CGC": "Arg",
    "CGA": "Arg", "CGG": "Arg",

    # AU
    "AUU": "Ile", "AUC": "Ile",
    "AUA": "Ile", "AUG": "Met",

    # AC
    "ACU": "Thr", "ACC": "Thr",
    "ACA": "Thr", "ACG": "Thr",

    # AA
    "AAU": "Asn", "AAC": "Asn",
    "AAA": "Lys", "AAG": "Lys",

    # AG
    "AGU": "Ser", "AGC": "Ser",
    "AGA": "Arg", "AGG": "Arg",

    # GU
    "GUU": "Val", "GUC": "Val",
    "GUA": "Val", "GUG": "Val",

    # GC
    "GCU": "Ala", "GCC": "Ala",
    "GCA": "Ala", "GCG": "Ala",

    # GA
    "GAU": "Asp", "GAC": "Asp",
    "GAA": "Glu", "GAG": "Glu",

    # GG
    "GGU": "Gly", "GGC": "Gly",
    "GGA": "Gly", "GGG": "Gly"
}

amino_acids = []
start_position = None

for i in range(len(mrna) - 2):
    if mrna[i:i + 3] == "AUG":
        start_position = i
        break

if start_position is None:
    print("Error: No start codon found.")

else:
    stop_codon = None

    for i in range(start_position, len(mrna) - 2, 3):
        codon = mrna[i:i + 3]

        if codon not in genetic_code:
            print("Unknown codon:", codon)
            amino_acids = []
            break

        if genetic_code[codon] == "Stop":
            stop_codon = codon
            break

        amino_acids.append(genetic_code[codon])

    if (len(mrna) - start_position) % 3 != 0 and stop_codon is None:
        print("Warning: Incomplete codon at the end.")

    print("mRNA sequence provided:", mrna)
    print("Start codon position:", start_position)
    print("Amino acid sequence:", "-".join(amino_acids))

    if stop_codon:
        print("Stop codon:", stop_codon)
    else:
        print("Warning: No stop codon found.")

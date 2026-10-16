from collections import Counter

DNA_CODON_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTC':'V', 'GTG':'V', 'GTT':'V', 'GTA':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'
}

def read_fasta(file_path):
    sequence_parts = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('>'):
                continue
            sequence_parts.append(line.strip().upper())
    return ''.join(sequence_parts)

def analyze_codons(file_path):
    genome = read_fasta(file_path)
    print(f"Total genome length: {len(genome):,} base pairs\n")
    
    codons = []
    for i in range(0, len(genome) - 2, 3):
        codon = genome[i:i+3]
        if all(base in 'ACGT' for base in codon):
            codons.append(codon)
            
    codon_counts = Counter(codons)
    total_codons = sum(codon_counts.values())
    
    top_10 = codon_counts.most_common(10)
    
    print("=" * 50)
    # table header
    print(f"{'Rank':<6}{'Codon':<8}{'RNA':<8}{'Amino Acid':<14}{'Frequency':<10}")
    print("=" * 50)
    
    for rank, (codon, count) in enumerate(top_10, 1):
        # Convert dna to rna by replacing T with U
        rna_codon = codon.replace('T', 'U')
        amino_acid = DNA_CODON_TABLE.get(codon, 'Unknown')
        
        print(f"{rank:<6}{codon:<8}{rna_codon:<8}{amino_acid:<14}{count:<10,}")
    print("=" * 50)

if __name__ == "__main__":
    file_name = "./GCF_001399695.1_ASM139969v1_cds_from_genomic.fna" 
    try:
        analyze_codons(file_name)
    except FileNotFoundError:
        print(f"Error: Could not find the file '{file_name}'")
import random

def generate_synthetic_genome():
 
    random.seed(42)
    bases = ['A', 'C', 'G', 'T']
    
    genome_list = [random.choice(bases) for _ in range(300)]    
    te1 = "ATCGATCGATCG"  # TE Alpha
    te2 = "ATCGAGCTAGCTA" # TE Beta
    te3 = "CCGATCGATCGG"  # TE Gamma
    te4 = "GCTAGCTAGCTA"  # TE Delta 
    
    # TE1 from index 50 to 62
    # TE2 from index 58 to 70 has an overlap
    # TE3 from index 200 to 211
    
    for i, base in enumerate(te1):
        genome_list[50 + i] = base

    for i, base in enumerate(te2):
        genome_list[58 + i] = base
        
    for i, base in enumerate(te3):
        genome_list[200 + i] = base

    for i, base in enumerate(te4):
        genome_list[250 + i] = base
        
    genome = "".join(genome_list)
    
    te_library = {
        "TE_Alpha": te1,
        "TE_Beta": te2,
        "TE_Gamma": te3,
        "TE_Delta": te4

    }
    
    return genome, te_library

def detect_transposons(genome, te_library):

    detections = []
    
    for te_name, te_seq in te_library.items():
        start_idx = 0
        # search for all occurrences of the TE sequence in the genome
        while True:
            pos = genome.find(te_seq, start_idx)
            if pos == -1:
                break
            end_idx = pos + len(te_seq)
            detections.append({
                "transposon": te_name,
                "sequence": te_seq,
                "start": pos,
                "end": end_idx
            })
            start_idx = pos + 1
            
    return detections

if __name__ == "__main__":
    print("Transposable Element detection simulation\n")
    
    genome, te_library = generate_synthetic_genome()
    print(f"Generated synthetic genome of length: {len(genome)} bases\n")
    
    print("Reference Transposons:")
    for name, seq in te_library.items():
        print(f"  - {name}: {seq} (Length: {len(seq)})")
    
    print("-" * 10)
    detected_elements = detect_transposons(genome, te_library)
    
    print(f"{'Transposon':<12}{'Start':<8}{'End':<8}{'Status':<15}")
    print("-" * 10)
    
    for item in detected_elements:
        print(f"{item['transposon']:<12}{item['start']:<8}{item['end']:<8}{'Detected':<15}")
        
    print("-" * 10)
    print("Overlap Analysis:")
    for i in range(len(detected_elements)):
        for j in range(i + 1, len(detected_elements)):
            t1 = detected_elements[i]
            t2 = detected_elements[j]
            
            overlap_start = max(t1["start"], t2["start"])
            overlap_end = min(t1["end"], t2["end"])
            
            if overlap_start < overlap_end:
                print(f"  -> Intersection detected between {t1['transposon']} "
                      f"[{t1['start']}-{t1['end']}] and {t2['transposon']} "
                      f"[{t2['start']}-{t2['end']}]. Overlapping span: indices {overlap_start} to {overlap_end}.")
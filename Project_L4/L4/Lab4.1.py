GENETIC_CODE = {
    'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
    'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
    'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'Stop', 'UAG': 'Stop',
    'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'Stop', 'UGG': 'Trp',

    'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
    'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
    'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',

    'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met (Start)',
    'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
    'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',

    'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
    'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
    'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
}

def dna_to_rna(dna_sequence: str) -> str:
    rna_sequence = dna_sequence.upper().replace('T', 'U')
    print(f"-> RNA Sequence (S_arn): {rna_sequence}")
    return rna_sequence


def rna_to_amino_acids(rna_sequence: str, genetic_code: dict) -> str:
    start_index = rna_sequence.find('AUG')

    if start_index == -1:
        return "Error: No 'AUG' start codon found in the RNA sequence."

    coding_rna = rna_sequence[start_index:]

    if len(coding_rna) % 3 != 0:
        coding_rna = coding_rna[:-(len(coding_rna) % 3)]

    amino_acid_sequence = []

    for i in range(0, len(coding_rna), 3):
        codon = coding_rna[i:i + 3]
        amino_acid = genetic_code.get(codon)

        if amino_acid is None:
            amino_acid_sequence.append(f"??? ({codon})")
            print(f"Warning: Found invalid codon {codon}. Stopping translation.")
            break

        if 'Stop' in amino_acid:
            print(f"-> Stop codon ({codon}) encountered. Translation terminated.")
            break

        amino_acid_sequence.append(amino_acid)

    return " - ".join(amino_acid_sequence)


def dna_to_protein(dna_sequence: str) -> dict:
    print(f"Initial DNA Sequence (S_adn): {dna_sequence}")

    rna_seq = dna_to_rna(dna_sequence)
    protein_seq = rna_to_amino_acids(rna_seq, GENETIC_CODE)

    return {
        "DNA_Sequence": dna_sequence,
        "RNA_Sequence": rna_seq,
        "Amino_Acid_Sequence": protein_seq
    }

example_dna = "CGAATGGGTATTCTAAGTTAGGTAAGAT"

print("GENETIC CODE TRANSLATOR")

results = dna_to_protein(example_dna)

print("\n FINAL RESULTS")
print(f"DNA Sequence:            {results['DNA_Sequence']}")
print(f"RNA Sequence:            {results['RNA_Sequence']}")
print(f"Amino Acid Sequence:     {results['Amino_Acid_Sequence']}")
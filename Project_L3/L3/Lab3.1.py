import math

def calculate_tm_basic(dna_sequence: str) -> float:
    dna_sequence = dna_sequence.upper()
    a = dna_sequence.count('A')
    t = dna_sequence.count('T')
    g = dna_sequence.count('G')
    c = dna_sequence.count('C')

    tm = 4 * (g + c) + 2 * (a + t)
    return tm


def calculate_tm_advanced(dna_sequence: str, na_conc: float = 0.05) -> float:
    dna_sequence = dna_sequence.upper()
    length = len(dna_sequence)
    g = dna_sequence.count('G')
    c = dna_sequence.count('C')
    gc_percent = (g + c) / length * 100

    tm = 81.5 + 16.6 * math.log10(na_conc) + 0.41 * gc_percent - 600 / length
    return tm


if __name__ == "__main__":
    dna = input("Enter DNA sequence (A, T, G, C): ").strip().upper()

    if not all(base in "ATGC" for base in dna):
        print("Error: DNA sequence must only contain A, T, G, or C.")
    else:
        tm_basic = calculate_tm_basic(dna)
        tm_advanced = calculate_tm_advanced(dna)

        print("\nResults:")
        print(f"DNA Sequence: {dna}")
        print(f"Length: {len(dna)} bases")
        print(f"GC%: {round((dna.count('G') + dna.count('C')) / len(dna) * 100, 2)}%")
        print(f"Basic Tm Formula: {tm_basic:.2f} °C")
        print(f"Advanced Tm Formula: {tm_advanced:.2f} °C")

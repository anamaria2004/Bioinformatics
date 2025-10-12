# DNA sequence
S = "ATTGTCCCAATCTGTTG"
bases = ["A", "C", "G", "T"]

# dinucleotides
dinucl = []
for nr1 in bases:
    for nr2 in bases:
        dinucl.append(nr1 + nr2)

#trinucleotides
trinucl = []
for nr1 in bases:
    for nr2 in bases:
        for nr3 in bases:
            trinucl.append(nr1 + nr2 + nr3)

#count overlapping occurrences
def overlapping_count(seq, sub):
    count = 0
    for i in range(len(seq) - len(sub) + 1):
        if seq[i:i+len(sub)] == sub:
            count += 1
    return count

#Brute-force engine for any k-mer
def brute_force_engine(S, leng, combos):
    total = len(S) - leng + 1
    print(f"{'Motif':<5} {'Count':<7} {'Frequency (%)':<15}")
    print("-" * 30)
    for c in combos:
        count = overlapping_count(S, c)
        percent = (count / total) * 100
        print(f"{c:<5} {count:<7} {percent:>10.2f}")
    print()

print("=== Dinucleotides (k=2) ===")
brute_force_engine(S, 2, dinucl)

print("=== Trinucleotides (k=3) ===")
brute_force_engine(S, 3, trinucl)

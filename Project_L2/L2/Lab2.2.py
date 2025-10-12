S = "ABAACVF"
alphabet = {"A", "C", "T", "G"}

#find existing dinucleotides
dinucl = set()
for i in range(len(S) - 1):
    pair = S[i:i+2]
    if all(base in alphabet for base in pair):
        dinucl.add(pair)

#find existing trinucleotides
trinucl = set()
for i in range(len(S) - 2):
    triple = S[i:i+3]
    if all(base in alphabet for base in triple):
        trinucl.add(triple)

print("=== Dinucleotides found ===")
print(sorted(dinucl))

print("\n=== Trinucleotides found ===")
print(sorted(trinucl))

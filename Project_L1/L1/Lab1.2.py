S = "ATTTCGCCGATA"

alphabet = []
for char in S:
    if char not in alphabet:
        alphabet.append(char)

freq = {}
for char in alphabet:
    count = 0
    for s in S:
        if s == char:
            count += 1
    freq[char] = count / len(S)

print("Alphabet:", alphabet)
print("Relative Frequencies:")
for char in alphabet:
    print(f"{char}: {freq[char]:.3f}")

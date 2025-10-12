S = "ATTTCGCCGATA"
alphabet=[]
for char in S:
    if char not in alphabet:
        alphabet.append(char)

print(alphabet)
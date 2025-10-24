import os
from Bio import SeqIO
from collections import Counter
import pandas as pd
import matplotlib.pyplot as plt
import textwrap

# --- Helper functions ---

def read_fasta_sequence(fasta_file):
    """Read FASTA file and concatenate all sequences."""
    records = list(SeqIO.parse(fasta_file, "fasta"))
    seq = "".join(str(r.seq).upper().replace("U", "T") for r in records)
    return seq

def count_codons(seq):
    """Count codons (triplets of A, T, G, C)."""
    codons = Counter()
    seq = seq.replace("\n", "").replace(" ", "")
    L = (len(seq) // 3) * 3
    for i in range(0, L, 3):
        codon = seq[i:i+3]
        if len(codon) == 3 and all(x in "ACGT" for x in codon):
            codons[codon] += 1
    return codons

genetic_code = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
    'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
    'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
    'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
    'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
    'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
    'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
    'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
    'GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}

def codon_to_aa_counts(codon_counts):
    """Convert codon frequencies to amino acid frequencies."""
    aa_counts = Counter()
    for codon, count in codon_counts.items():
        aa = genetic_code.get(codon, "X")
        aa_counts[aa] += count
    return aa_counts

def plot_top_codons(codon_counts, title, top_n=10):
    """Plot top codons."""
    top_items = codon_counts.most_common(top_n)
    codons, counts = zip(*top_items)
    plt.figure(figsize=(10,5))
    plt.bar(codons, counts, color='skyblue')
    plt.title(title)
    plt.xlabel("Codon")
    plt.ylabel("Frequency")
    for i, c in enumerate(counts):
        plt.text(i, c + max(counts)*0.01, str(c), ha='center', fontsize=9)
    plt.tight_layout()
    plt.show()

# --- Main program ---

print("🧬 Codon Frequency Analyzer")
print("---------------------------")

covid_file = input("Enter path to COVID-19 FASTA file: ").strip('"')
influenza_file = input("Enter path to Influenza FASTA file: ").strip('"')

if not os.path.exists(covid_file) or not os.path.exists(influenza_file):
    print("❌ One or both files not found. Please check the paths.")
    exit()

# Read and analyze
covid_seq = read_fasta_sequence(covid_file)
influenza_seq = read_fasta_sequence(influenza_file)

covid_codons = count_codons(covid_seq)
influenza_codons = count_codons(influenza_seq)

# --- Plot results ---
plot_top_codons(covid_codons, "Top 10 Codons – COVID-19 Genome")
plot_top_codons(influenza_codons, "Top 10 Codons – Influenza Genome")

# --- Compare common codons ---
covid_top10 = [c for c, _ in covid_codons.most_common(10)]
influenza_top10 = [c for c, _ in influenza_codons.most_common(10)]
common_codons = list(set(covid_top10) & set(influenza_top10))

print("\n🧩 Common codons in Top 10 of both genomes:")
print(common_codons if common_codons else "None found.")

# --- Amino acid frequencies ---
covid_aa = codon_to_aa_counts(covid_codons)
influenza_aa = codon_to_aa_counts(influenza_codons)

print("\n🔬 Top 3 amino acids – COVID-19:")
for aa, count in covid_aa.most_common(3):
    print(f"  {aa}: {count}")

print("\n🧬 Top 3 amino acids – Influenza:")
for aa, count in influenza_aa.most_common(3):
    print(f"  {aa}: {count}")

# --- AI Nutrition Prompt ---
covid_top3 = [aa for aa, _ in covid_aa.most_common(3)]
influenza_top3 = [aa for aa, _ in influenza_aa.most_common(3)]

prompt = f"""
You are a nutrition expert. Based on codon usage analysis:
- COVID-19's most common amino acids: {', '.join(covid_top3)}
- Influenza's most common amino acids: {', '.join(influenza_top3)}

Suggest 10 foods that are LOW in these amino acids.
For each, briefly explain why.
Include a mix of plant-based and common dietary items.
"""

print("\n💡 Suggested AI Prompt for further analysis:\n")
print(textwrap.fill(prompt.strip(), width=100))

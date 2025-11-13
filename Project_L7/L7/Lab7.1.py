import re
from collections import Counter, defaultdict
from Bio import Entrez, SeqIO
import matplotlib.pyplot as plt
import gradio as gr

Entrez.email = "your_email@example.com"

def download_sequence(accession):
    try:
        with Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text") as handle:
            record = SeqIO.read(handle, "fasta")
            return str(record.seq).upper()
    except Exception as e:
        return None

def detect_repeats(seq, min_len=3, max_len=6, min_repeats=2):
    repeats = defaultdict(Counter)
    for size in range(min_len, max_len + 1):
        for i in range(len(seq) - size + 1):
            motif = seq[i:i+size]
            repeats[size][motif] += 1

    filtered = {
        size: {motif: count for motif, count in counts.items() if count >= min_repeats}
        for size, counts in repeats.items()
    }
    return filtered

def plot_repeats(filtered):
    all_motifs = []
    for size, motifs in filtered.items():
        for motif, count in motifs.items():
            all_motifs.append((motif, count, size))

    if not all_motifs:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No significant repeats found.", ha='center', va='center')
        ax.axis("off")
        return fig

    all_motifs.sort(key=lambda x: -x[1])
    motifs, counts, sizes = zip(*all_motifs[:20])

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.barh(motifs, counts, color='royalblue')
    ax.set_xlabel("Repetition Count")
    ax.set_ylabel("Motif (3–6 bases)")
    ax.set_title("Top Repetitive DNA Motifs")
    plt.gca().invert_yaxis()
    plt.tight_layout()
    return fig

def analyze_dna_repeats(accession):
    seq = download_sequence(accession)
    if not seq:
        return "Error: Could not download sequence. Please check the accession ID.", None

    if len(seq) < 1000 or len(seq) > 3000:
        note = f"Warning: Sequence length = {len(seq)} (recommended 1000–3000 nt)"
    else:
        note = f"Sequence length = {len(seq)} nucleotides"

    filtered = detect_repeats(seq)
    fig = plot_repeats(filtered)

    summary = []
    for size, motifs in filtered.items():
        if motifs:
            top = sorted(motifs.items(), key=lambda x: -x[1])[:5]
            summary.append(f"Length {size} bp → {len(motifs)} unique repeats (Top: {top[0][0]} × {top[0][1]})")
    if not summary:
        summary_text = "No repetitive motifs (3–6 bp) found with at least 2 repetitions."
    else:
        summary_text = "\n".join(summary)

    return f"{note}\n\n{summary_text}", fig

iface = gr.Interface(
    fn=analyze_dna_repeats,
    inputs=gr.Textbox(label="Enter NCBI Accession ID (1000–3000 bp DNA sequence)", placeholder="e.g. NC_002695"),
    outputs=[gr.Textbox(label="Summary"), gr.Plot(label="Repetitive Motifs (3–6 bp)")],
    title="🧬 DNA Repetitive Sequence Detector",
    description="Fetches a DNA sequence from NCBI and detects repetitive motifs (3–6 bases) occurring at least twice."
)

iface.launch()
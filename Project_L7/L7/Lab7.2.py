import re
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np
from Bio import Entrez, SeqIO
import gradio as gr

Entrez.email = "your_email@example.com"

INFLUENZA_ACCESSIONS = [
    "CY121680", "CY121681", "CY121682", "CY121683", "CY121684",
    "CY121685", "CY121686", "CY121687", "CY121688", "CY121689"
]


def download_sequence(accession):
    with Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text") as handle:
        record = SeqIO.read(handle, "fasta")
        return str(record.seq).upper()


def find_repetitions(seq, min_len=3, max_len=6):
    repeats = []
    for size in range(min_len, max_len + 1):
        pattern = re.compile(rf'((\w{{{size}}}))\2{{1,}}')
        for match in pattern.finditer(seq):
            motif = match.group(2)
            repeats.append(motif)
    return repeats


def analyze_repeats(accessions):
    results = []
    for acc in accessions:
        seq = download_sequence(acc)
        repeats = find_repetitions(seq)
        if repeats:
            freq = Counter(repeats)
            top = freq.most_common(5)
            results.append({"accession": acc, "top": top})
        else:
            results.append({"accession": acc, "top": []})
    return results


def plot_repeats(results):
    fig, ax = plt.subplots(figsize=(12, 6))
    bar_width = 0.08
    positions = np.arange(len(results))

    for i, r in enumerate(results):
        motifs = [m for m, _ in r["top"]]
        counts = [c for _, c in r["top"]]
        ax.bar([i + j * bar_width for j in range(len(motifs))], counts, width=bar_width, label=r["accession"])
        for j, m in enumerate(motifs):
            ax.text(i + j * bar_width, counts[j] + 0.1, m, rotation=90, ha='center', fontsize=8)

    ax.set_xticks(positions)
    ax.set_xticklabels([r["accession"] for r in results], rotation=45)
    ax.set_ylabel("Repeat Count")
    ax.set_title("Most Frequent 3–6bp Repeats in Influenza Genomes")
    plt.tight_layout()
    return fig


def run_analysis():
    results = analyze_repeats(INFLUENZA_ACCESSIONS)
    fig = plot_repeats(results)
    summary = "\n".join([f"{r['accession']}: {r['top']}" for r in results])
    return summary, fig


iface = gr.Interface(
    fn=run_analysis,
    inputs=[],
    outputs=[gr.Textbox(label="Most Frequent Repeats Summary"), gr.Plot(label="Repeats Frequency Chart")],
    title="🦠 Influenza Genome Repeat Analyzer",
    description="Downloads 10 influenza genomes from NCBI and shows their most frequent DNA repeats."
)

iface.launch()
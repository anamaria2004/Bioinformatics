import random
import matplotlib.pyplot as plt
import numpy as np
from Bio import Entrez, SeqIO
import gradio as gr

Entrez.email = "your_email@example.com"


def simulate_gel(accession_or_seq="NC_012920.1", num_samples=10, min_len=100, max_len=3000):
    if all(c in "ATGCatgc\n " for c in accession_or_seq):
        dna_seq = accession_or_seq.upper().replace("\n", "").replace(" ", "")
    else:
        try:
            with Entrez.efetch(db="nucleotide", id=accession_or_seq, rettype="fasta", retmode="text") as handle:
                record = SeqIO.read(handle, "fasta")
                dna_seq = str(record.seq)
        except:
            return "Error fetching sequence. Check accession number.", None

    dna_seq = dna_seq[:3000]

    fragments = []
    for _ in range(num_samples):
        start = random.randint(0, len(dna_seq) - min_len)
        end = random.randint(start + min_len, min(len(dna_seq), start + max_len))
        fragments.append(dna_seq[start:end])
    fragment_lengths = [len(f) for f in fragments]

    log_lengths = np.log(fragment_lengths)
    max_distance = 10
    min_distance = 1
    distances = max_distance - (log_lengths - min(log_lengths)) / (max(log_lengths) - min(log_lengths)) * (
                max_distance - min_distance)

    fig, ax = plt.subplots(figsize=(6, 8))
    ax.set_title("Simulated DNA Gel Electrophoresis")
    ax.set_xlabel("Lane")
    ax.set_ylabel("Migration distance")

    for i, dist in enumerate(distances):
        ax.plot([i + 1] * 2, [0, dist], color="blue", linewidth=5, solid_capstyle="round")
        ax.text(i + 1, dist + 0.2, f"{fragment_lengths[i]} bp", ha="center")

    ax.set_xticks(range(1, num_samples + 1))
    ax.set_ylim(0, max_distance + 1)

    plt.tight_layout()
    return f"Fragment lengths: {fragment_lengths}", fig

iface = gr.Interface(
    fn=simulate_gel,
    inputs=[
        gr.Textbox(label="NCBI Accession Number or DNA Sequence", value="NC_012920.1"),
        gr.Slider(1, 20, value=10, step=1, label="Number of fragments"),
        gr.Slider(50, 5000, value=100, step=10, label="Minimum fragment length"),
        gr.Slider(100, 5000, value=3000, step=10, label="Maximum fragment length")
    ],
    outputs=[
        gr.Textbox(label="Fragment Information"),
        gr.Plot(label="Gel Simulation")
    ],
    title="DNA Gel Electrophoresis Simulator",
    description="Simulate migration of DNA fragments on a gel based on their lengths."
)

iface.launch()
import random
import re
import matplotlib.pyplot as plt
import numpy as np
from Bio import Entrez, SeqIO
import gradio as gr

Entrez.email = "your_email@example.com"


def download_sequence(accession):
    with Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text") as handle:
        record = SeqIO.read(handle, "fasta")
        return str(record.seq).upper()


def in_silico_digest(seq, enzyme_site="GAATTC"):
    fragments = re.split(enzyme_site, seq)
    lengths = [len(f) for f in fragments if len(f) > 0]
    return lengths


def simulate_gel_for_fragments(fragment_lengths, gel_max_dist=10, gel_min_dist=1):
    if len(fragment_lengths) == 0:
        return []
    log_lengths = np.log(fragment_lengths)
    min_log = np.min(log_lengths)
    max_log = np.max(log_lengths)
    if max_log == min_log:
        distances = np.full(len(fragment_lengths), (gel_max_dist + gel_min_dist) / 2)
    else:
        distances = gel_max_dist - (log_lengths - min_log) / (max_log - min_log) * (gel_max_dist - gel_min_dist)
    return distances


def plot_gel_columns(results):
    """
    Plot each genome as a separate column (like a gel lane),
    with dots representing fragment migration distances.
    """
    fig, ax = plt.subplots(figsize=(len(results) * 1.5, 8))

    for lane_idx, r in enumerate(results):
        for dist in r["distances"]:
            # x = lane index, y = migration distance
            ax.plot(lane_idx + 1, dist, 'o', color='blue', markersize=8)
        ax.text(lane_idx + 1, max(r["distances"]) + 0.5 if len(r["distances"]) > 0 else 0.5,
                r["accession"], rotation=90, ha='center', fontsize=8)

    ax.set_xticks(range(1, len(results) + 1))
    ax.set_xlabel("Genome (Lane)")
    ax.set_ylabel("Migration distance")
    ax.set_title("In-silico EcoRI Gel Simulation")
    ax.set_ylim(0, max([max(r["distances"]) if len(r["distances"]) > 0 else 0 for r in results]) + 2)
    plt.tight_layout()
    return fig


def process_genomes(accessions):
    results = []
    for acc in accessions:
        seq = download_sequence(acc)
        lengths = in_silico_digest(seq)
        distances = simulate_gel_for_fragments(lengths)
        results.append({"accession": acc, "lengths": lengths, "distances": distances})
    return results


def run_gui(accession_list_text):
    accessions = accession_list_text.strip().split()
    if len(accessions) < 10:
        return "Please provide 10 accessions separated by space.", None
    results = process_genomes(accessions[:10])

    # Determine which genome has most fragments
    max_frag = max(results, key=lambda r: len(r["lengths"]))
    summary = f"Genome with most fragments: {max_frag['accession']} ({len(max_frag['lengths'])} fragments)"

    fig = plot_gel_columns(results)
    return summary, fig


iface = gr.Interface(
    fn=run_gui,
    inputs=[gr.Textbox(label="Enter 10 accessions separated by space", value="")],
    outputs=[gr.Textbox(label="Summary"), gr.Plot(label="Gel Simulation")],
    title="Influenza Genome EcoRI Gel Simulation",
    description="Each genome is a separate column (lane). Simulates EcoRI digestion and fragment migration."
)

iface.launch()

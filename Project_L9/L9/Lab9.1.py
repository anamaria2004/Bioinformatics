import re
import gradio as gr
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

ENZYMES = {
    "EcoRI": {
        "seq": "GAATTC",
        "cut": 1
    },
    "BamHI": {
        "seq": "GGATCC",
        "cut": 1
    },
    "HindIII": {
        "seq": "AAGCTT",
        "cut": 1
    },
    "TaqI": {
        "seq": "TCGA",
        "cut": 1
    },
    "HaeIII": {
        "seq": "GGCC",
        "cut": 2
    }
}

def parse_fasta(raw_text):
    lines = raw_text.strip().split("\n")
    if lines[0].startswith(">"):
        lines = lines[1:]
    seq = "".join(lines).replace(" ", "").upper()
    seq = re.sub(r"[^ACGT]", "", seq)
    return seq

def digest_sequence(sequence, enzyme_name):
    enzyme = ENZYMES[enzyme_name]
    recog = enzyme["seq"]
    cut_offset = enzyme["cut"]

    cut_positions = []
    start = 0

    while True:
        idx = sequence.find(recog, start)
        if idx == -1:
            break
        cut_positions.append(idx + cut_offset)
        start = idx + 1

    cut_positions = sorted(list(set(cut_positions)))

    if not cut_positions:
        return [], [len(sequence)]

    fragments = []
    last = 0
    for cp in cut_positions:
        fragments.append(cp - last)
        last = cp
    fragments.append(len(sequence) - last)

    return cut_positions, fragments

def plot_gel(fragments_dict):
    fig, ax = plt.subplots(figsize=(4, 6))
    ax.set_facecolor("black")

    y_offset = 0
    for enzyme, frags in fragments_dict.items():
        sizes = sorted(frags, reverse=True)

        for size in sizes:
            y = np.log(size) * 20
            ax.hlines(y - y_offset, 0.2, 0.8, colors="white", linewidth=3)

        ax.text(0.9, -y_offset + 10, enzyme, color="white", fontsize=10)
        y_offset += 80

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("Simulated Electrophoresis Gel", color="white")
    plt.tight_layout()
    return fig

def process_sequence(raw_seq_text, selected_enzymes):

    seq = parse_fasta(raw_seq_text)
    length = len(seq)

    if length < 1000 or length > 3000:
        return (
            f"Sequence is {length} bp. Must be between 1000–3000 bp.",
            pd.DataFrame(),
            None
        )

    report = f"Sequence length: {length} bp\n"
    fragments_dict = {}
    rows = []

    for enzyme in selected_enzymes:
        cuts, frags = digest_sequence(seq, enzyme)
        fragments_dict[enzyme] = frags

        report += f"\n🔹 {enzyme}\n"
        report += f"Recognized sequence: {ENZYMES[enzyme]['seq']}\n"
        report += f"Cuts found: {len(cuts)}\n"
        report += f"Cut positions: {cuts if cuts else 'None'}\n"
        report += f"Fragment sizes: {frags}\n"

        for i, size in enumerate(frags):
            rows.append([enzyme, i+1, size])

    df = pd.DataFrame(rows, columns=["Enzyme", "Fragment #", "Length (bp)"])

    gel_fig = plot_gel(fragments_dict)

    return report, df, gel_fig

with gr.Blocks(theme=gr.themes.Soft(), title="Restriction Enzyme Simulator") as demo:

    gr.Markdown("""
    # Restriction Enzyme DNA Digestion Simulator  
    Paste any **1000–3000 bp DNA sequence** (FASTA-style allowed).  
    Select enzymes → Get digestion + gel simulation.  
    """)

    with gr.Row():
        seq_input = gr.Textbox(
            label="DNA Sequence (FASTA allowed)",
            placeholder="Paste FASTA sequence here...",
            lines=10
        )

        enzyme_select = gr.CheckboxGroup(
            choices=list(ENZYMES.keys()),
            value=["EcoRI", "BamHI"],
            label="Select Restriction Enzymes"
        )

    run_btn = gr.Button("Run Digestion", variant="primary")

    output_report = gr.Textbox(label="Results", lines=12)
    output_table = gr.Dataframe(label="Fragment Table")
    output_gel = gr.Plot(label="Simulated Gel")

    run_btn.click(
        process_sequence,
        inputs=[seq_input, enzyme_select],
        outputs=[output_report, output_table, output_gel]
    )

demo.launch()
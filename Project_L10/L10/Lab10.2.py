import gradio as gr
import numpy as np
import matplotlib.pyplot as plt
from Bio import Entrez, SeqIO
import io
import random

ENTREZ_EMAIL = "ncbi_user@example.com"

GENOMES_TO_ANALYZE = {
    "NC_045512.2": "COVID-19 | Wuhan-Hu-1",
    "OL440939.1": "COVID-19 | Delta",
    "OM005374.1": "COVID-19 | Omicron BA.1",
    "OP428801.1": "COVID-19 | Omicron BA.5",
    "OR593976.1": "COVID-19 | XBB.1.5",
    "ON696328.1": "COVID-19 | Kraken",
    "OR668962.1": "COVID-19 | EG.5",
    "OR873426.1": "COVID-19 | JN.1",
    "ON100742.1": "COVID-19 | Alpha",
    "OM922579.1": "COVID-19 | Gamma",

    "LC676834.1": "Influenza A | H1N1 (HA)",
    "LC733924.1": "Influenza A | H3N2 (HA)",
    "MW410023.1": "Influenza B | Victoria (HA)",
    "MW410029.1": "Influenza B | Yamagata (HA)",
    "MK038038.1": "Influenza A | H5N1 (HA)",
    "LC658804.1": "Influenza A | H1N1 (PB2)",
    "LC658805.1": "Influenza A | H1N1 (PB1)",
    "LC658806.1": "Influenza A | H1N1 (PA)",
    "MW410034.1": "Influenza B | Victoria (NA)",
    "MW410040.1": "Influenza B | Yamagata (NA)",
}

def fetch_genome_sequences(accession_ids):
    if not ENTREZ_EMAIL or ENTREZ_EMAIL == "your_email@example.com":
        return {"error": "Please set a valid NCBI email in the ENTREZ_EMAIL variable."}, {}

    Entrez.email = ENTREZ_EMAIL
    sequences = {}

    id_list = list(accession_ids.keys())
    chunks = [id_list[i:i + 5] for i in range(0, len(id_list), 5)]

    try:
        for chunk in chunks:
            handle = Entrez.efetch(db="nucleotide", id=chunk, rettype="fasta", retmode="text")
            fasta_data = handle.read()
            handle.close()

            for record in SeqIO.parse(io.StringIO(fasta_data), "fasta"):
                acc_id = record.id.split('|')[-1].split('.')[0]
                if acc_id in accession_ids:
                    sequences[accession_ids[acc_id]] = str(record.seq).upper().replace('N',
                                                                                       'A')

    except Exception as e:
        return {"error": f"An error occurred during NCBI retrieval: {e}"}, {}

    if len(sequences) != len(accession_ids):
        missing = set(accession_ids.values()) - set(sequences.keys())
        return {"warning": f"Could not retrieve all sequences. Missing: {list(missing)}"}, sequences

    return {}, sequences

def calculate_ods_coords(sequence, window_size=30):
    seq_len = len(sequence)
    if seq_len < window_size:
        return np.array([]), np.array([])

    ods_x = []
    ods_y = []

    for start in range(seq_len - window_size + 1):
        window = sequence[start:start + window_size]

        g_c_count = window.count('G') + window.count('C')
        a_t_count = window.count('A') + window.count('T')

        x_coord = g_c_count / window_size

        if g_c_count + a_t_count > 0:
            y_coord = (g_c_count / (g_c_count + a_t_count)) - (a_t_count / (g_c_count + a_t_count))
        else:
            y_coord = 0.0

        ods_x.append(x_coord)
        ods_y.append(y_coord)

    return np.array(ods_x), np.array(ods_y)


def calculate_center_of_weight(ods_x, ods_y):
    if len(ods_x) == 0:
        return 0, 0

    centroid_x = np.mean(ods_x)
    centroid_y = np.mean(ods_y)

    return centroid_x, centroid_y

def analyze_genomes_for_ods(genome_ids, window_size=30):

    status, sequences = fetch_genome_sequences(GENOMES_TO_ANALYZE)

    if status.get("error"):
        summary_text = status["error"]
        return summary_text, None, None

    all_ods_data = {}
    all_cw_data = {}

    for name, seq in sequences.items():
        ods_x, ods_y = calculate_ods_coords(seq, window_size)
        cw_x, cw_y = calculate_center_of_weight(ods_x, ods_y)

        all_ods_data[name] = (ods_x, ods_y)
        all_cw_data[name] = (cw_x, cw_y)

    fig1, ax1 = plt.subplots(figsize=(12, 8))

    ods_colors = {
        "COVID-19": 'blue',
        "Influenza": 'red'
    }

    for name, (ods_x, ods_y) in all_ods_data.items():
        virus_type = "COVID-19" if "COVID-19" in name else "Influenza"
        color = ods_colors[virus_type]

        ax1.scatter(ods_x, ods_y, alpha=0.1, s=5, color=color)

    ax1.set_title("Objective Digital Stains (ODS) for Influenza and COVID-19 Genomes")
    ax1.set_xlabel("Local GC Content (C+G Frequency)")
    ax1.set_ylabel("Local A/T vs C/G Bias Index")
    ax1.set_xlim(0, 1.0)
    ax1.legend(handles=[
        plt.Line2D([0], [0], marker='o', color='w', label='COVID-19 ODS', markerfacecolor='blue', markersize=10),
        plt.Line2D([0], [0], marker='o', color='w', label='Influenza ODS', markerfacecolor='red', markersize=10)
    ])
    ax1.grid(True, linestyle='--', alpha=0.5)
    plt.close(fig1)

    fig2, ax2 = plt.subplots(figsize=(10, 8))

    cw_x_all = []
    cw_y_all = []

    for name, (cw_x, cw_y) in all_cw_data.items():
        virus_type = "COVID-19" if "COVID-19" in name else "Influenza"
        color = ods_colors[virus_type]
        cw_x_all.append(cw_x)
        cw_y_all.append(cw_y)

        ax2.scatter(cw_x, cw_y, s=150, edgecolors=color, facecolors='none', linewidths=2)

        label_parts = name.split(' | ')
        label_text = f"{label_parts[0]}\n({label_parts[1]})"
        ax2.annotate(label_text, (cw_x, cw_y),
                     textcoords="offset points",
                     xytext=(5, 5),
                     ha='left',
                     fontsize=8,
                     color=color)

    ax2.set_title("Center of Weight (CW) of Objective Digital Stains")
    ax2.set_xlabel("CW X-Coordinate (Average GC Content)")
    ax2.set_ylabel("CW Y-Coordinate (Average A/T vs C/G Bias)")
    ax2.grid(True, linestyle=':', alpha=0.7)

    if cw_x_all and cw_y_all:
        margin_x = (max(cw_x_all) - min(cw_x_all)) * 0.1
        margin_y = (max(cw_y_all) - min(cw_y_all)) * 0.1
        ax2.set_xlim(min(cw_x_all) - margin_x, max(cw_x_all) + margin_x * 2.5)
        ax2.set_ylim(min(cw_y_all) - margin_y, max(cw_y_all) + margin_y * 1.5)

    ax2.legend(handles=[
        plt.Line2D([0], [0], marker='o', color='w', label='COVID-19 CW', markerfacecolor='none', markeredgecolor='blue',
                   markersize=10, linewidth=0),
        plt.Line2D([0], [0], marker='o', color='w', label='Influenza CW', markerfacecolor='none', markeredgecolor='red',
                   markersize=10, linewidth=0)
    ], loc='upper right')

    plt.close(fig2)

    # 5. Summary Text
    summary_text = (
        f"Successfully analyzed {len(sequences)} genomes (10 COVID-19, 10 Influenza) with a window size of {window_size}.\n"
        "The ODS plot shows the distribution of local genomic characteristics, and the Center of Weight (CW) plot provides a clear, labeled separation point for each genome/variant."
    )

    return summary_text, fig1, fig2


genome_id_list = "\n".join([f"{acc_id}: {name}" for acc_id, name in GENOMES_TO_ANALYZE.items()])

iface = gr.Interface(
    fn=analyze_genomes_for_ods,
    inputs=[
        gr.Textbox(label="Genome Accession IDs to Analyze (10 COVID-19, 10 Influenza)",
                   value=genome_id_list,
                   lines=20,
                   interactive=False),
        gr.Slider(minimum=5, maximum=100, step=5, value=30, label="Sliding Window Length (b) for ODS Calculation")
    ],
    outputs=[
        gr.Textbox(label="Analysis Status and Summary", lines=5),
        gr.Plot(label="Objective Digital Stains (ODS) - Local Frequency Pattern"),
        gr.Plot(label="Center of Weight (CW) of ODS - Labeled Classification")
    ],
    title=" Comparative Viral Genomics: Objective Digital Stain Analysis",
    description="Fetches 10 Influenza and 10 COVID-19 genomes from NCBI, plots their Objective Digital Stains (ODS), and maps the Centers of Weight (CW) for visual classification."
)

if __name__ == "__main__":
    iface.launch()

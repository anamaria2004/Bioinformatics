import gradio as gr
import numpy as np
import matplotlib.pyplot as plt
import requests
import io

def smith_waterman_score_native(seq1, seq2, match=3, mismatch=-3, gap=-2):
    rows = len(seq1) + 1
    cols = len(seq2) + 1

    previous_row = np.zeros(cols, dtype=int)
    current_row = np.zeros(cols, dtype=int)

    max_score = 0

    for i in range(1, rows):
        for j in range(1, cols):
            if seq1[i - 1] == seq2[j - 1]:
                diag = previous_row[j - 1] + match
            else:
                diag = previous_row[j - 1] + mismatch

            up = previous_row[j] + gap
            left = current_row[j - 1] + gap

            val = max(0, diag, up, left)
            current_row[j] = val

            if val > max_score:
                max_score = val

        previous_row[:] = current_row[:]

    return max_score

def fetch_ncbi_sequence(accession_id):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession_id}&rettype=fasta&retmode=text"
    try:
        response = requests.get(url)
        response.raise_for_status()
        data = response.text

        lines = data.strip().split('\n')
        if not lines:
            return None, "Error: Empty response"

        header = lines[0]
        sequence = "".join(lines[1:]).replace("\n", "").upper()
        return sequence, f"Successfully fetched {accession_id}\nHeader: {header[:50]}..."
    except Exception as e:
        return None, str(e)


def chunk_sequence(seq, size):
    return [seq[i:i + size] for i in range(0, len(seq), size)]


def run_windowed_alignment(seq1_id, seq2_id, window_size, progress=gr.Progress()):
    progress(0, desc="Fetching COVID-19 Genome...")
    s1, status1 = fetch_ncbi_sequence(seq1_id)
    if not s1: return status1, None

    progress(0.1, desc="Fetching Influenza Genome...")
    s2, status2 = fetch_ncbi_sequence(seq2_id)
    if not s2: return status2, None

    w_size = int(window_size)

    chunks1 = chunk_sequence(s1, w_size)
    chunks2 = chunk_sequence(s2, w_size)

    n_rows = len(chunks2)
    n_cols = len(chunks1)

    similarity_matrix = np.zeros((n_rows, n_cols))

    total_steps = n_rows * n_cols
    step_count = 0

    desc_str = f"Aligning {n_rows}x{n_cols} regions..."

    for i in range(n_rows):
        for j in range(n_cols):
            score = smith_waterman_score_native(chunks2[i], chunks1[j])

            similarity_matrix[i][j] = score

            step_count += 1
            if step_count % 50 == 0:
                progress((step_count / total_steps) * 0.9 + 0.1, desc=desc_str)

    fig, ax = plt.subplots(figsize=(12, 6))

    im = ax.imshow(similarity_matrix, cmap='inferno', aspect='auto', interpolation='nearest')

    ax.set_title(f"Genome Similarity Landscape\n(Window Size: {w_size}bp)")
    ax.set_xlabel(f"{seq1_id} (COVID-19) - Chunks")
    ax.set_ylabel(f"{seq2_id} (Influenza) - Chunks")

    cbar = plt.colorbar(im)
    cbar.set_label("Local Alignment Score (SW)")

    plt.tight_layout()

    result_text = (
        f"Analysis Complete.\n"
        f"Genome 1 ({seq1_id}): {len(s1)} bp -> {n_cols} regions\n"
        f"Genome 2 ({seq2_id}): {len(s2)} bp -> {n_rows} regions\n"
        f"Algorithm: Native Smith-Waterman on {w_size}bp windows.\n"
        f"Note: Dark regions indicate low similarity. Bright spots indicate conserved regions or high local alignment."
    )

    return result_text, fig


with gr.Blocks(title="Viral Genome Aligner") as app:
    gr.Markdown("## viral-align: Big Genome Local Alignment Tool")
    gr.Markdown(
        """
        This tool implements **Local Alignment (Smith-Waterman)** without external alignment libraries.
        To handle large files (like whole genomes) without crashing memory, it uses an **in-between layer** that breaks the genomes into 'windows' and aligns them grid-by-grid.
        """
    )

    with gr.Row():
        with gr.Column(scale=1):
            gr.Markdown("### Input Parameters")

            id1 = gr.Textbox(label="Genome 1 ID (NCBI)", value="NC_045512.2")
            id2 = gr.Textbox(label="Genome 2 ID (NCBI)", value="NC_026433.1")

            window_slider = gr.Slider(
                minimum=50,
                maximum=500,
                value=200,
                step=50,
                label="Window Size (bp)",
                info="Smaller windows = higher resolution but slower."
            )

            btn_align = gr.Button("Download & Align", variant="primary")

        with gr.Column(scale=2):
            output_text = gr.Textbox(label="Status / Logs", lines=5)
            output_plot = gr.Plot(label="Similarity Heatmap")

    btn_align.click(
        run_windowed_alignment,
        inputs=[id1, id2, window_slider],
        outputs=[output_text, output_plot]
    )

if __name__ == "__main__":
    app.launch()
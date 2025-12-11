import gradio as gr
import numpy as np
import matplotlib.pyplot as plt
import requests
import math

def smith_waterman_score_native(seq1, seq2, match=3, mismatch=-3, gap=-2):
    rows = len(seq1) + 1
    cols = len(seq2) + 1

    previous_row = [0] * cols
    current_row = [0] * cols
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

        previous_row = list(current_row)

    return max_score

def fetch_ncbi_sequence(accession_id):
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession_id}&rettype=fasta&retmode=text"
    try:
        response = requests.get(url)
        response.raise_for_status()
        lines = response.text.strip().split('\n')
        if not lines: return None, "Empty"
        return "".join(lines[1:]).replace("\n", "").upper(), "Success"
    except Exception as e:
        return None, str(e)


def chunk_sequence(seq, size):
    return [seq[i:i + size] for i in range(0, len(seq), size)]


def run_analysis(id1, id2, w_size, match_score, progress=gr.Progress()):
    progress(0, desc="Fetching...")
    s1, _ = fetch_ncbi_sequence(id1)
    s2, _ = fetch_ncbi_sequence(id2)

    if not s1 or not s2:
        return "Error fetching data", 0, 0, 0, None, None

    w_size = int(w_size)
    chunks1 = chunk_sequence(s1, w_size)
    chunks2 = chunk_sequence(s2, w_size)

    n_rows, n_cols = len(chunks2), len(chunks1)
    sim_matrix = np.zeros((n_rows, n_cols))

    total = n_rows * n_cols
    cnt = 0
    for i in range(n_rows):
        for j in range(n_cols):
            score = smith_waterman_score_native(chunks2[i], chunks1[j], match=match_score)
            sim_matrix[i][j] = score
            cnt += 1
            if cnt % 100 == 0:
                progress((cnt / total) * 0.8, desc="Aligning Windows...")

    total_score = np.sum(sim_matrix)
    min_len = min(len(s1), len(s2))
    max_theoretical = min_len * match_score
    s_norm = (total_score / max_theoretical) * 100 if max_theoretical > 0 else 0

    max_window_score = w_size * match_score
    thresh = max_window_score * 0.30
    high_conf = np.sum(sim_matrix > thresh)
    c_hc = (high_conf / total) * 100 if total > 0 else 0

    if total_score > 0:
        probs = sim_matrix.flatten() / total_score
        probs = probs[probs > 0]
        entropy = -np.sum(probs * np.log2(probs))
    else:
        entropy = 0

    fig1, ax1 = plt.subplots(figsize=(8, 5))
    ax1.imshow(sim_matrix, cmap='magma', aspect='auto')
    ax1.set_title("Similarity Heatmap")

    fig2, ax2 = plt.subplots(figsize=(8, 4))
    ax2.hist(sim_matrix.flatten(), bins=30, color='gray')
    ax2.axvline(thresh, color='red', linestyle='--')
    ax2.set_title("Score Distribution vs Threshold")

    status = f"Aligned {len(s1)}bp vs {len(s2)}bp in {n_cols}x{n_rows} windows."

    return status, round(s_norm, 3), round(c_hc, 3), round(entropy, 3), fig1, fig2

with gr.Blocks(title="Math-Based Genome Aligner") as app:
    gr.Markdown("## Genome Alignment Scoring")

    with gr.Row():
        # Inputs
        with gr.Column(scale=1):
            with gr.Group():
                inp_id1 = gr.Textbox(value="NC_045512.2", label="Sequence 1 (NCBI)")
                inp_id2 = gr.Textbox(value="NC_026433.1", label="Sequence 2 (NCBI)")
                inp_win = gr.Slider(50, 500, value=200, step=50, label="Window Size")
                inp_match = gr.Number(value=3, label="Match Score")
                btn = gr.Button("Calculate Scores", variant="primary")

            out_status = gr.Textbox(label="Status", interactive=False)

        with gr.Column(scale=2):
            gr.Markdown("### Mathematical Scoring Results")

            with gr.Row(variant="panel"):
                with gr.Column(scale=3):
                    gr.Markdown(
                        r"""
                        **1. Global Normalized Similarity**
                        $$S_{norm} = \frac{\sum S_{i,j}}{L_{min} \times Match} \times 100$$
                        *Measures total alignment mass relative to perfect identity.*
                        """
                    )
                with gr.Column(scale=1, min_width=100):
                    res_snorm = gr.Number(label="Result (%)")

            with gr.Row(variant="panel"):
                with gr.Column(scale=3):
                    gr.Markdown(
                        r"""
                        **2. High-Confidence Coverage**
                        $$C_{hc} = \frac{\text{Count}(S_{i,j} > 0.3 \times S_{max})}{N_{windows}} \times 100$$
                        *Percentage of regions with significant local similarity.*
                        """
                    )
                with gr.Column(scale=1, min_width=100):
                    res_chc = gr.Number(label="Result (%)")

            with gr.Row(variant="panel"):
                with gr.Column(scale=3):
                    gr.Markdown(
                        r"""
                        **3. Alignment Entropy**
                        $$E_{div} = - \sum p_i \log_2(p_i)$$
                        *Measures signal spread. Low = Focused homology, High = Random noise.*
                        """
                    )
                with gr.Column(scale=1, min_width=100):
                    res_entropy = gr.Number(label="Result")

    with gr.Row():
        out_plot1 = gr.Plot(label="Heatmap")
        out_plot2 = gr.Plot(label="Distribution")

    btn.click(
        run_analysis,
        inputs=[inp_id1, inp_id2, inp_win, inp_match],
        outputs=[out_status, res_snorm, res_chc, res_entropy, out_plot1, out_plot2]
    )

if __name__ == "__main__":
    app.launch()
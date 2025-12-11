import gradio as gr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors


def needleman_wunsch(seq1, seq2, match_score, mismatch_score, gap_penalty):
    n, m = len(seq1), len(seq2)

    score_matrix = np.zeros((m + 1, n + 1))

    for i in range(m + 1):
        score_matrix[i][0] = i * gap_penalty
    for j in range(n + 1):
        score_matrix[0][j] = j * gap_penalty

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq1[j - 1] == seq2[i - 1]:
                diagonal = score_matrix[i - 1][j - 1] + match_score
            else:
                diagonal = score_matrix[i - 1][j - 1] + mismatch_score

            up = score_matrix[i - 1][j] + gap_penalty
            left = score_matrix[i][j - 1] + gap_penalty

            score_matrix[i][j] = max(diagonal, up, left)

    align1, align2 = "", ""
    i, j = m, n
    path_coords = [(i, j)]
    matches = 0

    while i > 0 or j > 0:
        current_score = score_matrix[i][j]

        if i > 0 and j > 0:
            if seq1[j - 1] == seq2[i - 1]:
                score_diag = score_matrix[i - 1][j - 1] + match_score
                is_match = True
            else:
                score_diag = score_matrix[i - 1][j - 1] + mismatch_score
                is_match = False

            if current_score == score_diag:
                align1 += seq1[j - 1]
                align2 += seq2[i - 1]
                if is_match:
                    matches += 1
                i -= 1
                j -= 1
                path_coords.append((i, j))
                continue

        if i > 0 and current_score == score_matrix[i - 1][j] + gap_penalty:
            align1 += "-"
            align2 += seq2[i - 1]
            i -= 1
            path_coords.append((i, j))
            continue

        if j > 0 and current_score == score_matrix[i][j - 1] + gap_penalty:
            align1 += seq1[j - 1]
            align2 += "-"
            j -= 1
            path_coords.append((i, j))
            continue

    return score_matrix, align1[::-1], align2[::-1], path_coords, matches


def process_alignment(s1, s2, gap, match, mismatch):
    matrix, a1, a2, path, matches = needleman_wunsch(s1, s2, match, mismatch, gap)

    length = len(a1)
    similarity = int((matches / length) * 100) if length > 0 else 0

    bars = ""
    for k in range(length):
        if a1[k] == a2[k] and a1[k] != "-":
            bars += "|"
        else:
            bars += " "

    text_result = (
        f"{a1}\n"
        f"{bars}\n"
        f"{a2}\n\n"
        f"Matches = {matches}\n"
        f"Length = {length}\n"
        f"Similarity = {similarity} %"
    )

    fig_heatmap, ax1 = plt.subplots(figsize=(6, 5))
    im = ax1.imshow(matrix, cmap='magma', aspect='auto')
    ax1.set_title("Graphic representation of alignment matrix")
    ax1.set_xlabel("Sequence 1")
    ax1.set_ylabel("Sequence 2")
    fig_heatmap.colorbar(im)
    plt.tight_layout()

    fig_grid, ax2 = plt.subplots(figsize=(6, 5))
    rows, cols = matrix.shape

    grid_display = np.ones((rows, cols, 3))
    grid_display[:, :] = [1, 1, 0.85]

    for r, c in path:
        grid_display[r, c] = [0.8, 0.2, 0.2]  # Red

    ax2.imshow(grid_display)

    ax2.set_xticks(np.arange(-0.5, cols, 1), minor=True)
    ax2.set_yticks(np.arange(-0.5, rows, 1), minor=True)
    ax2.grid(which='minor', color='black', linestyle='-', linewidth=0.5)
    ax2.tick_params(which="minor", size=0)

    if len(s1) < 30 and len(s2) < 30:
        ax2.set_xticks(range(cols))
        ax2.set_xticklabels(['-'] + list(s1))
        ax2.set_yticks(range(rows))
        ax2.set_yticklabels(['-'] + list(s2))

    ax2.set_title("Traceback path")
    plt.tight_layout()

    return text_result, fig_heatmap, fig_grid

with gr.Blocks(title="DNA Alignment App") as app:
    gr.Markdown("## DNA Sequence Alignment (Needleman-Wunsch)")

    with gr.Row():
        with gr.Column(scale=1):
            s1_input = gr.Textbox(
                label="Sequence 1",
                value="ACCGTGAAGCCAATAC"
            )
            s2_input = gr.Textbox(
                label="Sequence 2",
                value="AGCGTGCAGCCAATAC"
            )

            with gr.Group():
                gr.Markdown("### Parameters")
                gap_input = gr.Number(label="Gap", value=0)
                match_input = gr.Number(label="Match", value=1)
                mismatch_input = gr.Number(label="Mismatch", value=-1)

            btn = gr.Button("Align", variant="primary")

        with gr.Column(scale=2):
            with gr.Row():
                heatmap_out = gr.Plot(label="Alignment Matrix")
                grid_out = gr.Plot(label="Traceback Path")

            text_out = gr.Textbox(
                label="Alignment Result",
                lines=6,
                elem_classes="output-text",
                show_copy_button=True,
                elem_id="mono_font"
            )

    app.css = "#mono_font textarea { font-family: 'Courier New', monospace; }"

    btn.click(
        process_alignment,
        inputs=[s1_input, s2_input, gap_input, match_input, mismatch_input],
        outputs=[text_out, heatmap_out, grid_out]
    )

if __name__ == "__main__":
    app.launch()
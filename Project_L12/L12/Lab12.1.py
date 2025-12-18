import numpy as np
import pandas as pd
import gradio as gr

bases = ['A', 'C', 'G', 'T']
counts = {
    1: [3, 2, 1, 4], 2: [6, 2, 1, 1], 3: [1, 1, 7, 1],
    4: [0, 0, 10, 0], 5: [0, 0, 0, 10], 6: [6, 2, 1, 1],
    7: [7, 1, 1, 1], 8: [2, 1, 5, 2], 9: [1, 2, 1, 6]
}

count_df = pd.DataFrame(counts, index=bases)
freq_df = count_df / 10
null_model = 0.25

epsilon = 0.0001
ll_matrix_df = np.log(freq_df.replace(0, epsilon) / null_model).round(3)


def run_analysis(sequence):
    s_clean = sequence.upper().strip().replace(" ", "")
    window_size = 9
    results = []

    for i in range(len(s_clean) - window_size + 1):
        window = s_clean[i:i + window_size]
        score = 0
        for pos, char in enumerate(window):
            if char in bases:
                score += ll_matrix_df.loc[char, pos + 1]
            else:
                score -= 10
        results.append({"Index": i, "Window": window, "Score": round(score, 3)})

    res_df = pd.DataFrame(results).sort_values(by="Score", ascending=False)

    c_out = count_df.reset_index().rename(columns={'index': 'Base'})
    f_out = freq_df.reset_index().rename(columns={'index': 'Base'})
    l_out = ll_matrix_df.reset_index().rename(columns={'index': 'Base'})

    return c_out, f_out, l_out, res_df


with gr.Blocks(title="Exon-Intron Motif Finder") as demo:
    gr.Markdown("# DNA Splice Site Recognition")
    gr.Markdown("Implementation of Count, Frequency, and Log-Likelihood analysis for Motif Finding.")

    with gr.Row():
        input_seq = gr.Textbox(
            label="Input Sequence S",
            value="CAGGTTGGAAACGTAATCAGCGATTACGCATGACGTAA",
            lines=2
        )

    btn = gr.Button("Generate Matrices & Analyze", variant="primary")

    with gr.Tabs():
        with gr.TabItem("1. Count Matrix"):
            gr.Markdown("Absolute occurrences of each nucleotide per position.")
            out_count = gr.DataFrame()

        with gr.TabItem("2. Frequency/Weight Matrix"):
            gr.Markdown("Relative frequencies (P(N)) normalized by n=10.")
            out_freq = gr.DataFrame()

        with gr.TabItem("3. Log-Likelihoods"):
            gr.Markdown("Values calculated using ln(P(N)/0.25).")
            out_llm = gr.DataFrame()

        with gr.TabItem("4. Sequence Analysis"):
            gr.Markdown("Sliding window scores for the provided sequence S.")
            out_res = gr.DataFrame()

    btn.click(
        fn=run_analysis,
        inputs=input_seq,
        outputs=[out_count, out_freq, out_llm, out_res]
    )

if __name__ == "__main__":
    demo.launch()
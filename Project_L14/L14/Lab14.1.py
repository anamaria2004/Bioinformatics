import gradio as gr
import numpy as np
import pandas as pd
import math

S1_POSITIVE = "ATCGATTCGATATCATACACGTAT"
S2_NEGATIVE = "CTCGACTAGTATGAAGTCCACGCTTG"
BASES = ['A', 'C', 'G', 'T']

def count_transitions(sequence, bases):
    counts = {b1: {b2: 1 for b2 in bases} for b1 in bases}
    for i in range(len(sequence) - 1):
        curr = sequence[i]
        nxt = sequence[i + 1]
        if curr in bases and nxt in bases:
            counts[curr][nxt] += 1
    return counts


def calculate_probabilities(counts, bases):
    probs = {b1: {b2: 0.0 for b2 in bases} for b1 in bases}
    for b1 in bases:
        row_sum = sum(counts[b1].values())
        for b2 in bases:
            probs[b1][b2] = counts[b1][b2] / row_sum
    return probs


def create_log_likelihood_matrix(p_plus, p_minus, bases):
    ll_matrix = {b1: {b2: 0.0 for b2 in bases} for b1 in bases}
    for b1 in bases:
        for b2 in bases:
            val_plus = p_plus[b1][b2]
            val_neg = p_minus[b1][b2]
            ll_matrix[b1][b2] = math.log2(val_plus / val_neg)
    return ll_matrix


def format_df(data_dict):
    df = pd.DataFrame.from_dict(data_dict, orient='index')
    df = df[BASES].reindex(BASES)
    return df.round(4)

counts_pos_dict = count_transitions(S1_POSITIVE, BASES)
counts_neg_dict = count_transitions(S2_NEGATIVE, BASES)

probs_pos_dict = calculate_probabilities(counts_pos_dict, BASES)
probs_neg_dict = calculate_probabilities(counts_neg_dict, BASES)

ll_dict = create_log_likelihood_matrix(probs_pos_dict, probs_neg_dict, BASES)

df_counts_pos = format_df(counts_pos_dict)
df_counts_neg = format_df(counts_neg_dict)
df_probs_pos = format_df(probs_pos_dict)
df_probs_neg = format_df(probs_neg_dict)
df_ll = format_df(ll_dict)

def analyze_sequence(sequence):
    sequence = sequence.upper().strip()

    if not all(c in BASES for c in sequence):
        return "Error: Invalid chars", 0.0

    score = 0.0
    for i in range(len(sequence) - 1):
        c1 = sequence[i]
        c2 = sequence[i + 1]
        score += df_ll.loc[c1, c2]

    decision = "✅ CpG Island (Positive)" if score > 0 else "❌ Non-Island (Negative)"
    return decision, round(score, 4)


with gr.Blocks(title="CpG Island Detector", theme=gr.themes.Soft()) as demo:
    gr.Markdown("# CpG Island Detector: Full Pipeline Analysis")
    gr.Markdown(f"**Training Data:**\n* Positive Model (S1): `{S1_POSITIVE}`\n* Negative Model (S2): `{S2_NEGATIVE}`")

    with gr.Row():
        with gr.Column():
            inp_seq = gr.Textbox(label="Input DNA Sequence", value="CAGGTTGGAAACGTAA")
            btn_run = gr.Button("Run Analysis", variant="primary")
        with gr.Column():
            out_res = gr.Textbox(label="Final Decision")
            out_score = gr.Number(label="Total Log-Likelihood Score")

    gr.Markdown("### Intermediate Calculation Steps")

    with gr.Tabs():
        with gr.TabItem("1. Count Matrices"):
            gr.Markdown("Raw transition counts (with Laplace smoothing +1 applied).")
            with gr.Row():
                gr.Dataframe(value=df_counts_pos, label="CpG+ Counts (S1)", interactive=False)
                gr.Dataframe(value=df_counts_neg, label="CpG- Counts (S2)", interactive=False)

        with gr.TabItem("2. Probability Matrices"):
            gr.Markdown("Normalized probabilities (Row Sum = 1).")
            with gr.Row():
                gr.Dataframe(value=df_probs_pos, label="CpG+ Probabilities", interactive=False)
                gr.Dataframe(value=df_probs_neg, label="CpG- Probabilities", interactive=False)

        with gr.TabItem("3. Log-Likelihood Matrix"):
            gr.Dataframe(value=df_ll, label="Log-Likelihood Matrix (Bits)", interactive=False)

    btn_run.click(fn=analyze_sequence, inputs=inp_seq, outputs=[out_res, out_score])

if __name__ == "__main__":
    demo.launch()
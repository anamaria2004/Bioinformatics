import gradio as gr
import pandas as pd
import numpy as np
import math
import re
import matplotlib.pyplot as plt
from collections import defaultdict

EMINESCU_TEXT = """
A fost odată ca-n poveşti,
A fost ca niciodată,
Din rude mari împărăteşti,
O prea frumoasă fată.
Şi era una la părinţi
Şi mândră-n toate cele,
Cum e Fecioara între sfinţi
Şi luna între stele.
Cobori în jos, luceafăr blând,
Alunecând pe-o rază,
Pătrunde-n casă şi în gând
Şi viaţa-mi luminează!
"""

STANESCU_TEXT = """
Leoaica tânără, iubirea
mi-a sărit în faţă.
Mă pândise-n încordare
mai demult.
Colţii albi mi i-a înfipt în faţă,
m-a muşcat leoaica, azi, de faţă.
Şi deodată-n jurul meu, natura
se făcu un cerc, de-a-dura,
când mai larg, când mai aproape,
ca o strângere de ape.
Şi privirea-n sus ţâşni,
curcubeu tăiat în două,
şi auzul o-ntâlni
tocmai lângă ciocârlii.
"""

SUSPECT_TEXT = """
A fost odată ca-n poveşti o leoaică tânără.
Cobori în jos luceafăr blând şi colţii albi mi i-a înfipt în faţă.
Din rude mari împărăteşti privirea-n sus ţâşni
şi era una la părinţi un cerc de-a-dura.
"""

def clean_and_tokenize(text):
    text = re.sub(r'[^\w\s]', '', text)
    tokens = text.lower().split()
    return tokens

def train_markov_model(text):
    tokens = clean_and_tokenize(text)
    transitions = defaultdict(lambda: defaultdict(int))

    for i in range(len(tokens) - 1):
        curr_w, next_w = tokens[i], tokens[i + 1]
        transitions[curr_w][next_w] += 1

    model = defaultdict(dict)
    vocab = set(tokens)

    for curr_w in transitions:
        total_counts = sum(transitions[curr_w].values())
        for next_w in vocab:
            count = transitions[curr_w].get(next_w, 0)
            prob = (count + 1) / (total_counts + len(vocab))
            model[curr_w][next_w] = prob

    return model, vocab

def calculate_llr(model_a, model_b, vocab_a, vocab_b, sequence):
    tokens = clean_and_tokenize(sequence)
    results = []
    all_vocab_len = len(vocab_a.union(vocab_b))
    default_prob = 1 / (all_vocab_len * 100)

    for i in range(len(tokens) - 1):
        curr_w = tokens[i]
        next_w = tokens[i + 1]

        if curr_w in model_a:
            p_a = model_a[curr_w].get(next_w, default_prob)
        else:
            p_a = default_prob

        if curr_w in model_b:
            p_b = model_b[curr_w].get(next_w, default_prob)
        else:
            p_b = default_prob

        try:
            llr = math.log(p_a / p_b)
        except ValueError:
            llr = 0.0

        if abs(llr) < 0.1:
            llr = 0.0

        results.append((curr_w, next_w, llr))

    return results


def analyze_plagiarism(ref_eminescu, ref_stanescu, suspect_text):
    model_e, vocab_e = train_markov_model(ref_eminescu)
    model_s, vocab_s = train_markov_model(ref_stanescu)

    transitions = calculate_llr(model_e, model_s, vocab_e, vocab_s, suspect_text)

    html_output = "<div style='line-height: 2.5; font-size: 18px;'>"

    eminescu_count = 0
    stanescu_count = 0
    neutral_count = 0

    y_values = []
    x_positions = []
    current_char_pos = 0
    raw_scores = []

    for i, (w1, w2, score) in enumerate(transitions):
        current_char_pos += len(w1) + 1
        x_positions.append(current_char_pos)
        raw_scores.append(score)

        if score > 0.5:
            color = "#d1fae5"
            border = "2px solid #10b981"
            tooltip = f"Eminescu (Score: {score:.2f})"
            eminescu_count += 1
        elif score < -0.5:
            color = "#fee2e2"
            border = "2px solid #ef4444"
            tooltip = f"Stanescu (Score: {score:.2f})"
            stanescu_count += 1
        else:
            color = "#f3f4f6"
            border = "1px dashed #9ca3af"
            tooltip = "Neutral / Common"
            neutral_count += 1

        span = f"""
        <span style="background-color: {color}; color: black; border-bottom: {border}; padding: 2px 5px; margin: 0 2px; border-radius: 4px;" title="{tooltip}">
            {w1} {w2}
        </span>
        """
        html_output += span

    html_output += "</div>"

    window_size = 4
    scores_series = pd.Series(raw_scores)
    smooth_scores = scores_series.rolling(window=window_size, min_periods=1, center=True).sum()

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.plot(x_positions, smooth_scores, color='#1f77b4', linewidth=2, label='Style Trajectory')

    ax.axhline(0, color='#0f172a', linewidth=1.5)

    ax.set_title("Eminescu (above 0) vs Stănescu (below 0)", fontsize=12, fontweight='bold')
    ax.set_xlabel("Text position (chars)")
    ax.set_ylabel("Log-likelihood ratio (Window Sum)")
    ax.grid(True, linestyle='--', alpha=0.5)

    ax.fill_between(x_positions, smooth_scores, 0, where=(smooth_scores > 0), color='green', alpha=0.1)
    ax.fill_between(x_positions, smooth_scores, 0, where=(smooth_scores < 0), color='red', alpha=0.1)

    plt.tight_layout()

    total = len(transitions)
    if total == 0: total = 1
    stats = {
        "Eminescu %": round(eminescu_count / total * 100, 1),
        "Stanescu %": round(stanescu_count / total * 100, 1),
        "Neutral %": round(neutral_count / total * 100, 1)
    }

    return html_output, pd.DataFrame([stats]), fig

custom_css = """
.header {text-align: center; background: #374151; color: white; padding: 20px; border-radius: 10px; margin-bottom: 20px;}
.header h1 {margin: 0; font-size: 2em;}
.legend {margin-top: 10px; font-weight: bold;}
.eminescu-tag {color: #059669;} 
.stanescu-tag {color: #dc2626;}
"""

with gr.Blocks(theme=gr.themes.Soft(), css=custom_css, title="Courtroom AI") as demo:
    with gr.Column():
        gr.HTML("""
        <div class="header">
            <h1>The AI Courtroom: Plagiarism Detector</h1>
            <p>Distinguishing styles using Markov Chain Probabilities</p>
        </div>
        """)

        with gr.Row():
            with gr.Column():
                gr.Markdown("### Model A: Mihai Eminescu (Reference)")
                inp_eminescu = gr.Textbox(value=EMINESCU_TEXT, lines=8, label="Training Corpus A")

            with gr.Column():
                gr.Markdown("###Model B: Nichita Stanescu (Reference)")
                inp_stanescu = gr.Textbox(value=STANESCU_TEXT, lines=8, label="Training Corpus B")

        gr.Markdown("---")

        with gr.Row():
            with gr.Column(scale=1):
                gr.Markdown("### The Accused Text (Mihai)")
                gr.Markdown("_This text is a generated mix of both authors to test the detector._")
                inp_suspect = gr.Textbox(value=SUSPECT_TEXT, lines=5, label="Suspect Text")
                btn = gr.Button("Analyze Style (Judge)", variant="primary")

            with gr.Column(scale=1):
                gr.Markdown("### Forensic Graph Analysis")
                out_plot = gr.Plot(label="Sliding Window Log-Likelihood")
                out_stats = gr.Dataframe(label="Overall Stats")

        gr.Markdown("### Detailed Text Map")
        gr.HTML("""
        <div class="legend">
            Legend: 
            <span class="eminescu-tag">■ Green = Eminescu Style</span> | 
            <span class="stanescu-tag">■ Red = Stanescu Style</span> | 
            <span style="color: grey">■ Grey = Neutral/Unknown</span>
        </div>
        """)
        out_html = gr.HTML(label="Highlighted Text")

    btn.click(fn=analyze_plagiarism,
              inputs=[inp_eminescu, inp_stanescu, inp_suspect],
              outputs=[out_html, out_stats, out_plot])

if __name__ == "__main__":
    demo.launch()
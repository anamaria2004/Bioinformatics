import gradio as gr
import random
import json
import pandas as pd
import numpy as np


def analyze_transitions(length=50):
    bases = ['A', 'C', 'G', 'T']
    sequence = "".join(random.choices(bases, k=length))

    counts = {b1: {b2: 0 for b2 in bases} for b1 in bases}

    for i in range(len(sequence) - 1):
        current_base = sequence[i]
        next_base = sequence[i + 1]
        counts[current_base][next_base] += 1

    transition_matrix = {b1: {b2: 0.0 for b2 in bases} for b1 in bases}

    for base in bases:
        total_transitions = sum(counts[base].values())
        if total_transitions > 0:
            for next_base in bases:
                prob = counts[base][next_base] / total_transitions
                transition_matrix[base][next_base] = round(prob, 3)
        else:
            pass

    output_data = {
        "sequence_metadata": {
            "length": len(sequence),
            "sequence": sequence
        },
        "transition_matrix": transition_matrix
    }

    json_filename = "transition_matrix.json"
    with open(json_filename, "w") as f:
        json.dump(output_data, f, indent=4)

    df = pd.DataFrame.from_dict(transition_matrix, orient='index')

    return sequence, df, json_filename

with gr.Blocks(title="DNA Transition Probability Calculator") as demo:
    gr.Markdown("# DNA Transition Matrix Calculator")
    gr.Markdown("Generates a random DNA sequence and calculates the probability of each nucleotide following another.")

    with gr.Row():
        generate_btn = gr.Button("Generate Random Sequence & Calculate", variant="primary")

    with gr.Row():
        with gr.Column():
            gr.Markdown("### 1. Generated Sequence")
            seq_output = gr.Textbox(label="Random Sequence (50bp)", lines=3)

            gr.Markdown("### 3. Download Result")
            file_output = gr.File(label="Download JSON")

        with gr.Column():
            gr.Markdown("### 2. Transition Matrix (Probabilities)")
            matrix_output = gr.DataFrame(label="Transition Probabilities")

    generate_btn.click(
        fn=analyze_transitions,
        inputs=[],
        outputs=[seq_output, matrix_output, file_output]
    )

if __name__ == "__main__":
    demo.launch()
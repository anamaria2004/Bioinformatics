import gradio as gr
import json
import random
import pandas as pd


def load_and_generate(file_obj, length=50, start_word=""):
    if file_obj is None:
        return "Please upload a JSON file first."

    try:
        with open(file_obj.name, 'r') as f:
            data = json.load(f)

        mapping = data.get("mapping", {})
        matrix = data.get("transition_matrix", {})
        symbol_to_word = {v: k for k, v in mapping.items()}

    except Exception as e:
        return f"Error reading JSON: {str(e)}"

    current_symbol = None

    if start_word:
        clean_start = start_word.lower().strip()
        if clean_start in mapping:
            current_symbol = mapping[clean_start]
        else:
            return f"Error: The word '{start_word}' is not in the source vocabulary."
    else:
        current_symbol = random.choice(list(matrix.keys()))

    generated_symbols = [current_symbol]

    for _ in range(int(length) - 1):
        if current_symbol not in matrix:
            break

        transitions = matrix[current_symbol]

        next_symbols = list(transitions.keys())
        probabilities = list(transitions.values())

        if not next_symbols:
            break

        next_symbol = random.choices(next_symbols, weights=probabilities, k=1)[0]

        generated_symbols.append(next_symbol)
        current_symbol = next_symbol

    output_words = [symbol_to_word.get(sym, "???") for sym in generated_symbols]
    output_text = " ".join(output_words)

    return output_text

with gr.Blocks(title="Text Synthesizer (Markov Chain)") as demo:
    gr.Markdown("# Text Synthesizer from Transition Matrix")
    gr.Markdown("Upload the `word_transitions.json` file created in the previous step to generate new text.")

    with gr.Row():
        file_input = gr.File(label="Upload 'word_transitions.json'", file_types=[".json"])

    with gr.Row():
        with gr.Column():
            length_slider = gr.Slider(minimum=5, maximum=100, value=20, step=1, label="Number of Words to Generate")
            start_box = gr.Textbox(label="Starting Word (Optional)", placeholder="Leave empty for random start")
            gen_btn = gr.Button("Synthesize Text", variant="primary")

        with gr.Column():
            output_box = gr.Textbox(label="Generated Text", lines=5)

    gen_btn.click(
        fn=load_and_generate,
        inputs=[file_input, length_slider, start_box],
        outputs=[output_box]
    )

if __name__ == "__main__":
    demo.launch()
import gradio as gr
import random
import json
import pandas as pd
import string
from collections import defaultdict

SENTENCES = [
    "The quick brown fox jumps over the lazy dog.",
    "Data science is an interdisciplinary field.",
    "Artificial intelligence is transforming the world.",
    "Machine learning models require good data.",
    "Python is a powerful programming language.",
    "Gradio makes it easy to build user interfaces.",
    "Transition matrices are used in Markov chains.",
    "Probability is the measure of the likelihood that an event will occur.",
    "Deep learning is a subset of machine learning.",
    "Natural language processing helps computers understand human language."
]


def generate_text(target_length=300):
    text = ""
    while len(text) < target_length:
        text += random.choice(SENTENCES) + " "
    return text.strip()[:target_length + 50]


def analyze_word_transitions(custom_text=None):
    if not custom_text:
        text = generate_text()
    else:
        text = custom_text

    translator = str.maketrans('', '', string.punctuation)
    tokens = [word.translate(translator).lower() for word in text.split()]
    tokens = [t for t in tokens if t]

    unique_words = sorted(list(set(tokens)))

    word_to_symbol = {}
    symbol_to_word = {}

    start_ascii = 33
    for i, word in enumerate(unique_words):
        char_code = start_ascii + (i % 90)
        symbol = chr(char_code)
        if i >= 90:
            symbol += str(i // 90)

        word_to_symbol[word] = symbol
        symbol_to_word[symbol] = word

    transitions = defaultdict(lambda: defaultdict(int))

    for i in range(len(tokens) - 1):
        curr_word = tokens[i]
        next_word = tokens[i + 1]

        curr_sym = word_to_symbol[curr_word]
        next_sym = word_to_symbol[next_word]

        transitions[curr_sym][next_sym] += 1

    matrix_dict = {}
    all_symbols = list(symbol_to_word.keys())

    for sym_from in all_symbols:
        matrix_dict[sym_from] = {}
        total = sum(transitions[sym_from].values())

        for sym_to in all_symbols:
            if total > 0:
                count = transitions[sym_from].get(sym_to, 0)
                matrix_dict[sym_from][sym_to] = round(count / total, 3)
            else:
                matrix_dict[sym_from][sym_to] = 0.0

    output_data = {
        "text": text,
        "mapping": word_to_symbol,
        "transition_matrix": matrix_dict
    }

    filename = "word_transitions.json"
    with open(filename, 'w') as f:
        json.dump(output_data, f, indent=4)

    mapping_data = [{"Symbol": k, "Word": v} for k, v in symbol_to_word.items()]
    df_mapping = pd.DataFrame(mapping_data)

    df_matrix = pd.DataFrame.from_dict(matrix_dict, orient='index')
    df_matrix = df_matrix.fillna(0.0)

    return text, df_mapping, df_matrix, filename

with gr.Blocks(title="Word Transition Probability (Markov)") as demo:
    gr.Markdown("# Word Transition Matrix Calculator")
    gr.Markdown("Generates English text, maps words to ASCII symbols, and calculates transition probabilities.")

    with gr.Row():
        input_text = gr.Textbox(label="Input Text (Optional - leave empty for random)", lines=2)
        btn = gr.Button("Analyze Text", variant="primary")

    with gr.Row():
        text_display = gr.Textbox(label="Analyzed Text (~300 chars)", lines=3)

    with gr.Row():
        with gr.Column(scale=1):
            gr.Markdown("### Symbol Mapping")
            map_output = gr.DataFrame(label="Word -> Symbol")

        with gr.Column(scale=2):
            gr.Markdown("### Transition Matrix (Probabilities)")
            matrix_output = gr.DataFrame(label="P(Next Symbol | Current Symbol)")

    file_dl = gr.File(label="Download JSON Result")

    btn.click(
        fn=analyze_word_transitions,
        inputs=[input_text],
        outputs=[text_display, map_output, matrix_output, file_dl]
    )

if __name__ == "__main__":
    demo.launch()
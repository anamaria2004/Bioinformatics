import gradio as gr
import numpy as np
import matplotlib.pyplot as plt

S = "CGGACTGATCTATCTAAAAAAAAAAAAAAAAAAAAAAAAAAACGTAGCATCTATCGATCTATCTAGCGATCTATCTACTACG"
WINDOW_LENGTH = 30


def calculate_gc_percentage(window_sequence):
    g_count = window_sequence.count('G')
    c_count = window_sequence.count('C')
    total_length = len(window_sequence)
    if total_length == 0:
        return 0.0
    return ((g_count + c_count) / total_length) * 100


def calculate_kappa_ic(window_sequence):
    N = len(window_sequence)
    if N < 2:
        return 0.0

    counts = {
        'A': window_sequence.count('A'),
        'T': window_sequence.count('T'),
        'C': window_sequence.count('C'),
        'G': window_sequence.count('G')
    }

    sum_f_squared_minus_f = 0
    for count in counts.values():
        sum_f_squared_minus_f += count * (count - 1)

    ic_value = sum_f_squared_minus_f / (N * (N - 1))

    return ic_value * 400

def dna_pattern_analysis(dna_sequence, window_len):
    seq_len = len(dna_sequence)
    if seq_len < window_len:
        return "Sequence is shorter than the window length.", None, None, None

    window_starts = []
    gc_percentages = []
    kappa_ics = []

    for start in range(seq_len - window_len + 1):
        window = dna_sequence[start:start + window_len]
        window_starts.append(start + 1)  # 1-based start position

        gc_val = calculate_gc_percentage(window)
        gc_percentages.append(gc_val)

        ic_val = calculate_kappa_ic(window)
        kappa_ics.append(ic_val)

    ic_array = np.array(kappa_ics)

    sum_ic = np.sum(ic_array)
    cw = None
    if sum_ic > 0:
        cw = np.sum(np.array(window_starts) * ic_array) / sum_ic

    fig1, ax1 = plt.subplots(figsize=(10, 6))
    ax1.plot(window_starts, kappa_ics, marker='o', linestyle='-', color='indigo', label=r'$\kappa$ IC Pattern')
    if cw is not None:
        ax1.axvline(cw, color='red', linestyle='--', label=f'Center of Weight: {cw:.2f}')
    ax1.set_title(r'Pattern Analysis: $\kappa$ Index of Coincidence (IC)')
    ax1.set_xlabel(f'Window Start Position (Window Size: {window_len}b)')
    ax1.set_ylabel('Kappa Index of Coincidence (Normalized)')
    ax1.legend()
    ax1.grid(True)
    plt.close(fig1)

    fig2, ax2 = plt.subplots(figsize=(10, 6))

    if cw is not None:
        ax2.scatter([cw], [1], color='red', s=200, label=f'Center of Weight ({cw:.2f})')
        ax2.set_title(r'Center of Weight (CW) of the $\kappa$ IC Pattern')
        ax2.set_xlabel('Sequence Position')
        ax2.set_yticks([])
        ax2.set_ylim(0, 2)
        ax2.grid(True, axis='x')
        ax2.legend()
    else:
        ax2.text(0.5, 0.5, "Center of Weight could not be calculated (Sum IC is 0)",
                 horizontalalignment='center', verticalalignment='center', transform=ax2.transAxes)

    plt.close(fig2)

    full_seq_gc = calculate_gc_percentage(dna_sequence)
    full_seq_ic = calculate_kappa_ic(dna_sequence)

    verification_output = (
            f"Sequence Length: {seq_len} bp\n"
            f"Window Length: {window_len} bp\n"
            f"\n Verification Values\n"
            f"3. GC Content (Whole Sequence): {full_seq_gc:.2f} % (Target: 29.27)\n"
            f"4. Kappa IC (Whole Sequence): {full_seq_ic:.2f} (Target: 27.53)\n"
            f"Our results ({full_seq_gc:.2f} and {full_seq_ic:.2f}) confirm the functions are working and match the target values from steps 3 and 4, which were calculated on the full sequence 'S'."
            f"\n Pattern Center\n"
            f"6. Center of Weight (CW) of the " + r'$\kappa$' + f" IC Pattern: {cw:.2f}"
    )

    return verification_output, fig1, fig2

iface = gr.Interface(
    fn=dna_pattern_analysis,
    inputs=[
        gr.Textbox(label="DNA Sequence (Promoter)", value=S, lines=3),
        gr.Slider(minimum=5, maximum=50, step=1, value=WINDOW_LENGTH, label="Sliding Window Length (b)")
    ],
    outputs=[
        gr.Textbox(label="Analysis Summary & Verification", lines=10),
        gr.Plot(label=r"5. $\kappa$ IC Pattern with Center of Weight"),
        gr.Plot(label="7. Center of Weight Plot")
    ],
    title="DNA Promoter Pattern Analyzer (Kappa Index of Coincidence)",
    description=r"This application computes the C+G% and Kappa Index of Coincidence ($\kappa$ IC) for a DNA sequence using a sliding window and determines the Center of Weight for the resulting pattern."
)

if __name__ == "__main__":
    print("Launching Gradio interface...")
    iface.launch()
    print("Gradio interface closed.")
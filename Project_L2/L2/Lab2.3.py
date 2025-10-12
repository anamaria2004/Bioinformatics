import tkinter as tk
from tkinter import filedialog, messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

WINDOW_SIZE = 30
ALPHABET = ("A", "C", "G", "T")


def read_fasta_sequence(path):
    seq_parts = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                continue
            seq_parts.append(line.upper())
    return "".join(seq_parts)


def prefix_counts(sequence):
    n = len(sequence)
    prefixes = {b: np.zeros(n + 1, dtype=np.int32) for b in ALPHABET}

    for i, ch in enumerate(sequence):
        for b in ALPHABET:
            prefixes[b][i + 1] = prefixes[b][i] + (1 if ch == b else 0)

    return prefixes


def compute_window_frequencies(sequence, window_size=WINDOW_SIZE):
    n = len(sequence)
    if n < window_size:
        return [], {b: np.array([]) for b in ALPHABET}

    prefixes = prefix_counts(sequence)
    total_windows = n - window_size + 1
    positions = np.arange(total_windows, dtype=int)

    freqs = {}
    for b in ALPHABET:
        counts = prefixes[b][window_size:] - prefixes[b][: -window_size]
        freqs[b] = counts.astype(float) / window_size

    return positions, freqs


class FastaWindowApp:
    def __init__(self, master):
        self.master = master
        master.title("FASTA Sliding Window Frequency Viewer")

        self.frame_controls = tk.Frame(master)
        self.frame_controls.pack(fill=tk.X, padx=8, pady=6)

        self.btn_open = tk.Button(self.frame_controls, text="Select FASTA file", command=self.select_file)
        self.btn_open.pack(side=tk.LEFT)

        self.lbl_info = tk.Label(self.frame_controls, text="No file selected.")
        self.lbl_info.pack(side=tk.LEFT, padx=10)

        self.fig, self.ax = plt.subplots(figsize=(8, 4))
        self.ax.set_xlabel("Window start index")
        self.ax.set_ylabel("Relative frequency (0-1)")
        self.ax.set_title(f"Sliding window (size {WINDOW_SIZE}) nucleotide frequencies")

        self.canvas = FigureCanvasTkAgg(self.fig, master=master)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(fill=tk.BOTH, expand=True, padx=8, pady=6)

        # status
        self.status = tk.Label(master, text="", bd=1, relief=tk.SUNKEN, anchor=tk.W)
        self.status.pack(fill=tk.X)

    def select_file(self):
        file_path = filedialog.askopenfilename(
            title="Choose FASTA file",
            filetypes=[("FASTA files", "*.fasta *.faa *.fna"), ("All files", "*.*")]
        )
        if not file_path:
            return
        try:
            self.status.config(text="Reading file...")
            self.master.update_idletasks()

            seq = read_fasta_sequence(file_path)
            seq_len = len(seq)
            self.lbl_info.config(text=f"File: {file_path}  |  Sequence length: {seq_len}")

            if seq_len == 0:
                messagebox.showerror("Empty sequence", "No biological sequence (A/C/G/T) found in the file.")
                self.status.config(text="No sequence found.")
                return

            self.status.config(text="Computing sliding-window frequencies...")
            self.master.update_idletasks()
            positions, freqs = compute_window_frequencies(seq, WINDOW_SIZE)

            if len(positions) == 0:
                messagebox.showinfo("Short sequence", f"Sequence length < {WINDOW_SIZE}; no windows to analyze.")
                self.status.config(text="Sequence too short.")
                self.clear_plot()
                return

            self.ax.clear()
            self.ax.set_title(f"Sliding window (size {WINDOW_SIZE}) nucleotide frequencies")
            self.ax.set_xlabel("Window start index")
            self.ax.set_ylabel("Relative frequency (0-1)")
            for b in ALPHABET:
                self.ax.plot(positions, freqs[b], label=b)  # do NOT set explicit colors

            self.ax.legend(title="Nucleotide")
            self.ax.set_ylim(0, 1)
            self.ax.grid(True, linestyle=':', linewidth=0.5)
            self.canvas.draw()

            self.status.config(text=f"Done — {len(positions)} windows analyzed.")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to process file:\n{e}")
            self.status.config(text="Error during processing.")

    def clear_plot(self):
        self.ax.clear()
        self.canvas.draw()


def run_app():
    root = tk.Tk()
    app = FastaWindowApp(root)
    root.geometry("900x600")
    root.mainloop()


if __name__ == "__main__":
    run_app()

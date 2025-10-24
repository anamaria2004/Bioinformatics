import tkinter as tk
from tkinter import filedialog, messagebox
import matplotlib.pyplot as plt
from math import log10
from Bio import SeqIO
import os

WINDOW_SIZE = 9
NA_CONCENTRATION = 0.05

def tm_simple(seq):
    """Formula 1: Tm = 4(G + C) + 2(A + T)"""
    g = seq.count("G")
    c = seq.count("C")
    a = seq.count("A")
    t = seq.count("T")
    return 4 * (g + c) + 2 * (a + t)


def tm_advanced(seq, na_conc=NA_CONCENTRATION):
    """Formula 2: Tm = 81.5 + 16.6*log10([Na+]) + 0.41*(%GC) – 600/length"""
    length = len(seq)
    gc_percent = 100 * (seq.count("G") + seq.count("C")) / length
    return 81.5 + 16.6 * log10(na_conc) + 0.41 * gc_percent - (600 / length)


def sliding_window(seq, window_size):
    for i in range(len(seq) - window_size + 1):
        yield seq[i:i + window_size]


def analyze_sequence(fasta_path):
    record = next(SeqIO.parse(fasta_path, "fasta"))
    dna_seq = str(record.seq).upper().replace("N", "")
    tm1_values, tm2_values, positions = [], [], []

    for i, window in enumerate(sliding_window(dna_seq, WINDOW_SIZE)):
        tm1_values.append(tm_simple(window))
        tm2_values.append(tm_advanced(window))
        positions.append(i + 1)

    return positions, tm1_values, tm2_values


def plot_tm_chart(positions, tm1, tm2, title):
    plt.figure(figsize=(10, 5))
    plt.plot(positions, tm1, label="Formula 1: 4(G+C)+2(A+T)", color="blue")
    plt.plot(positions, tm2, label="Formula 2: Advanced (Na+)", color="red")
    plt.title(title)
    plt.xlabel("Window Start Position")
    plt.ylabel("Melting Temperature (°C)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()


class DNAApp:
    def __init__(self, master):
        self.master = master
        master.title("DNA Melting Temperature Analyzer")
        master.geometry("500x250")
        master.configure(bg="#f0f4f8")

        self.label = tk.Label(master, text="Upload a FASTA file to analyze:", font=("Arial", 12), bg="#f0f4f8")
        self.label.pack(pady=15)

        self.upload_button = tk.Button(master, text="📂 Upload FASTA", command=self.upload_file, font=("Arial", 11), bg="#0078D7", fg="white", padx=10, pady=5)
        self.upload_button.pack()

        self.file_label = tk.Label(master, text="", font=("Arial", 10), bg="#f0f4f8", fg="#333")
        self.file_label.pack(pady=10)

        self.analyze_button = tk.Button(master, text="🔬 Analyze Sequence", command=self.run_analysis, font=("Arial", 11), bg="#28a745", fg="white", padx=10, pady=5, state=tk.DISABLED)
        self.analyze_button.pack()

        self.fasta_path = None

    def upload_file(self):
        file_path = filedialog.askopenfilename(filetypes=[("FASTA Files", "*.fasta *.faa")])
        if file_path:
            self.fasta_path = file_path
            self.file_label.config(text=os.path.basename(file_path))
            self.analyze_button.config(state=tk.NORMAL)

    def run_analysis(self):
        if not self.fasta_path:
            messagebox.showwarning("No file", "Please upload a FASTA file first.")
            return

        try:
            positions, tm1, tm2 = analyze_sequence(self.fasta_path)
            title = f"Melting Temperature: {os.path.basename(self.fasta_path)} (window={WINDOW_SIZE})"
            plot_tm_chart(positions, tm1, tm2, title)
        except Exception as e:
            messagebox.showerror("Error", str(e))


if __name__ == "__main__":
    root = tk.Tk()
    app = DNAApp(root)
    root.mainloop()

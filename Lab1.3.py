import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext

def process_fasta(file_path):
    freq = {"A": 0, "C": 0, "G": 0, "T": 0}  # only DNA bases
    total = 0

    with open(file_path, "r") as f:
        header = f.readline()  # read the header line

        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            for char in line:
                if char in freq:
                    freq[char] += 1
                    total += 1

    if total == 0:
        raise ValueError("No valid DNA bases (A, C, G, T) found in sequence.")

    relative_freq = {k: v / total for k, v in freq.items()}
    return header.strip(), freq, relative_freq, total

def open_file():
    file_path = filedialog.askopenfilename(
        title="Select FASTA File",
        filetypes=(("FASTA files", "*.fasta *.faa *.fna"), ("All files", "*.*"))
    )

    if file_path:
        try:
            header, freq, rel_freq, total = process_fasta(file_path)

            result_text.delete(1.0, tk.END)
            result_text.insert(tk.END, f"Header: {header}\n")
            result_text.insert(tk.END, f"Sequence length (A,C,G,T only): {total:,}\n\n")
            result_text.insert(tk.END, "Alphabet and Relative Frequencies:\n")
            for char in ["A", "C", "G", "T"]:  # consistent order
                result_text.insert(tk.END, f"{char}: {rel_freq[char]:.4f}\n")

        except Exception as e:
            messagebox.showerror("Error", f"Failed to process file:\n{e}")

root = tk.Tk()
root.title("FASTA DNA Analyzer")
root.geometry("500x400")

open_button = tk.Button(root, text="Choose FASTA File", command=open_file)
open_button.pack(pady=10)

result_text = scrolledtext.ScrolledText(root, wrap=tk.WORD, width=60, height=20)
result_text.pack(padx=10, pady=10, fill=tk.BOTH, expand=True)

root.mainloop()

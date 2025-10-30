import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext
from Bio import SeqIO
import random

def upload_fasta():
    filepath = filedialog.askopenfilename(
        title="Select FASTA file",
        filetypes=[("FASTA files", "*.fasta *.fa *.txt"), ("All files", "*.*")]
    )
    if not filepath:
        return
    fasta_path.set(filepath)
    text_box.insert(tk.END, f"\n Loaded FASTA file: {filepath}\n")

def process_sequence():
    path = fasta_path.get()
    if not path:
        messagebox.showerror("Error", "Please upload a FASTA file first.")
        return

    try:
        record = next(SeqIO.parse(path, "fasta"))
    except Exception as e:
        messagebox.showerror("Error", f"Failed to parse FASTA: {e}")
        return

    seq = str(record.seq).upper()
    seq_len = len(seq)

    text_box.insert(tk.END, f"\n Original sequence length: {seq_len}\n")
    text_box.insert(tk.END, f"Original (first 300 nt): {seq[:300]}\n\n")

    if seq_len < 1000 or seq_len > 3000:
        text_box.insert(tk.END, f" Warning: Sequence length not in 1000–3000 range\n\n")

    samples = []
    for _ in range(2000):
        length = random.randint(100, 150)
        start = random.randint(0, seq_len - length)
        fragment = seq[start:start + length]
        samples.append((start, fragment))

    text_box.insert(tk.END, f" Generated {len(samples)} random samples (100–150 nt each)\n")
    text_box.insert(tk.END, "Example samples (with positions):\n")
    for i, (pos, frag) in enumerate(samples[:5]):
        text_box.insert(tk.END, f"  Sample {i+1}: start={pos}, seq={frag[:80]}...\n")
    text_box.insert(tk.END, "\n")

    text_box.insert(tk.END, " Rebuilding original sequence by sorting fragments...\n")
    samples.sort(key=lambda x: x[0])

    reconstructed = ["N"] * seq_len
    for start, frag in samples:
        reconstructed[start:start+len(frag)] = list(frag)
    reconstructed_seq = "".join(reconstructed)

    text_box.insert(tk.END, f" Reconstruction complete.\n")
    text_box.insert(tk.END, f"Reconstructed sequence length: {len(reconstructed_seq)}\n")
    text_box.insert(tk.END, f"Reconstructed (first 300 nt): {reconstructed_seq[:300]}\n\n")

    if reconstructed_seq == seq:
        text_box.insert(tk.END, " The reconstructed sequence is IDENTICAL to the original!\n", "success")
    else:
        match_count = sum(1 for a, b in zip(seq, reconstructed_seq) if a == b)
        similarity = match_count / len(seq) * 100
        text_box.insert(tk.END, f" The sequences differ.\n", "error")
        text_box.insert(tk.END, f"Similarity: {similarity:.2f}%\n")

    text_box.insert(tk.END, "\n" + "="*70 + "\n")

root = tk.Tk()
root.title("DNA Sequence Rebuilder GUI")
root.geometry("900x600")

fasta_path = tk.StringVar()

top_frame = tk.Frame(root)
top_frame.pack(pady=10)

upload_btn = tk.Button(top_frame, text="Upload FASTA File", command=upload_fasta, bg="#d0ebff")
upload_btn.pack(side=tk.LEFT, padx=5)

process_btn = tk.Button(top_frame, text="Process Sequence", command=process_sequence, bg="#d8f5a2")
process_btn.pack(side=tk.LEFT, padx=5)

text_box = scrolledtext.ScrolledText(root, wrap=tk.WORD, width=110, height=30, font=("Courier", 9))
text_box.pack(padx=10, pady=10)

text_box.tag_config("success", foreground="green", font=("Courier", 10, "bold"))
text_box.tag_config("error", foreground="red", font=("Courier", 10, "bold"))

root.mainloop()
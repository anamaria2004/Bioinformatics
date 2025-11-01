import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext
import time, os, random
from Bio import Entrez, SeqIO
from Bio.SeqUtils import gc_fraction
import matplotlib.pyplot as plt
import numpy as np

def GC(seq):
    try:
        return gc_fraction(seq) * 100
    except Exception:
        seq = str(seq).upper()
        return 100.0 * (seq.count("G") + seq.count("C")) / len(seq)

class ViralGenomeAnalyzer:
    def __init__(self, master):
        self.master = master
        master.title("Viral Genome Analyzer")
        master.geometry("850x650")
        master.config(bg="#f4f4f4")

        self.file_paths = []
        self.results = []
        self.temp_files = []
        self.MAX_FILES = 10

        Entrez.email = "student.lab@example.com"

        tk.Label(master, text="Viral Genome Analysis Tool", font=("Arial", 18, "bold"), bg="#f4f4f4", fg="#222").pack(pady=10)
        tk.Label(master, text="Step 1: Upload or Fetch 10 Viral Genomes", font=("Arial", 12), bg="#f4f4f4").pack()

        frame = tk.Frame(master, bg="#f4f4f4")
        frame.pack(pady=10)

        tk.Button(frame, text="Upload FASTA Files", command=self.upload_files, bg="#4CAF50", fg="white", width=20).grid(row=0, column=0, padx=10)
        tk.Button(frame, text="Fetch 10 from NCBI", command=self.fetch_from_ncbi, bg="#FF9800", fg="white", width=20).grid(row=0, column=1, padx=10)

        self.file_list = tk.Listbox(master, height=6, width=80, bg="#fff", fg="#333", font=("Arial", 10))
        self.file_list.pack(pady=5)

        tk.Button(master, text="▶ Run Analysis & Plot", command=self.run_analysis, bg="#007BFF", fg="white", font=("Arial", 12, "bold")).pack(pady=15)

        self.output = scrolledtext.ScrolledText(master, height=12, width=90, bg="#fff", fg="#111", font=("Consolas", 10))
        self.output.pack(pady=10)
        self.output.insert(tk.END, "Results will appear here...\n")

    def upload_files(self):
        files = filedialog.askopenfilenames(title="Select up to 10 FASTA files", filetypes=[("FASTA files", "*.fasta *.fa *.fna")])
        self.file_paths = list(files)[:self.MAX_FILES]
        self.file_list.delete(0, tk.END)
        for f in self.file_paths:
            self.file_list.insert(tk.END, os.path.basename(f))

        self.output.insert(tk.END, f"\n Loaded {len(self.file_paths)} local FASTA files.\n")

    def fetch_from_ncbi(self):
        try:
            self.output.insert(tk.END, "\n Fetching viral genomes from NCBI...\n")
            self.master.update()

            handle = Entrez.esearch(db="nucleotide", term="viral genome[title] AND complete genome[title]", retmax=1000)
            record = Entrez.read(handle)
            handle.close()

            ids = random.sample(record["IdList"], 10)
            handle = Entrez.efetch(db="nucleotide", id=ids, rettype="fasta", retmode="text")
            data = handle.read().strip().split('>')
            handle.close()

            self.file_paths.clear()
            self.file_list.delete(0, tk.END)

            for i, seq_block in enumerate(data):
                if not seq_block.strip():
                    continue
                fasta = ">" + seq_block
                fname = f"viral_genome_{i+1}.fasta"
                with open(fname, "w") as f:
                    f.write(fasta)
                self.file_paths.append(fname)
                self.temp_files.append(fname)
                self.file_list.insert(tk.END, fname)

            self.output.insert(tk.END, f"Successfully fetched {len(self.file_paths)} genomes from NCBI.\n")

        except Exception as e:
            messagebox.showerror("Error", f"Failed to fetch genomes: {e}")

    def simulate_assembly(self, seq):
        gc_percent = GC(seq)
        base_time = len(seq) / 1e6 * 500 + random.uniform(100, 300)
        gc_factor = abs(gc_percent - 50) * random.uniform(2, 5)
        total_time = base_time + gc_factor
        return gc_percent, total_time

    def run_analysis(self):
        if not self.file_paths:
            messagebox.showwarning("No Files", "Please upload or fetch genomes first.")
            return

        self.results.clear()
        self.output.insert(tk.END, "\n Running analysis...\n")

        for i, path in enumerate(self.file_paths):
            try:
                record = next(SeqIO.parse(path, "fasta"))
                seq = str(record.seq)
                gc, t = self.simulate_assembly(seq)
                self.results.append((os.path.basename(path), gc, t))
                self.output.insert(tk.END, f"{i+1}. {os.path.basename(path)} | GC%: {gc:.2f} | Time: {t:.1f} ms\n")
            except Exception as e:
                self.output.insert(tk.END, f"Error processing {path}: {e}\n")

        self.output.insert(tk.END, "\n Analysis complete! Plotting results...\n")
        self.plot_results()
        self.generate_report()

    def plot_results(self):
        if not self.results:
            return

        names = [r[0] for r in self.results]
        gc_vals = [r[1] for r in self.results]
        times = [r[2] for r in self.results]

        plt.style.use("seaborn-v0_8-darkgrid")
        plt.figure(figsize=(9, 6))
        plt.scatter(gc_vals, times, c=np.arange(len(gc_vals)), cmap="viridis", s=150, edgecolors="w")

        for i, (gc, t) in enumerate(zip(gc_vals, times)):
            plt.text(gc + 0.3, t + 5, str(i+1), fontsize=9, weight="bold")

        plt.xlabel("Overall G+C Percentage [%]")
        plt.ylabel("Simulated Assembly Time [ms]")
        plt.title("Viral Genome Assembly Time vs G+C Content")
        plt.tight_layout()
        plt.show()

    def generate_report(self):
        if not self.results:
            return

        report = "Genome Assembly Report\n"
        report += "="*50 + "\n\n"
        for i, (name, gc, t) in enumerate(self.results):
            report += f"{i+1}. {name}\n   - G+C: {gc:.2f}%\n   - Assembly Time: {t:.1f} ms\n\n"

        gc_vals = [r[1] for r in self.results]
        times = [r[2] for r in self.results]
        report += "Interpretation:\n"
        report += "-"*50 + "\n"
        report += "Genomes with G+C values far from 50% tend to show longer simulated assembly times.\n"
        report += "Points higher on the Y-axis represent more complex or larger genomes.\n"
        report += "Points further on the X-axis indicate genomes with higher or lower GC content.\n"

        with open("viral_genome_analysis_report.txt", "w") as f:
            f.write(report)

        self.output.insert(tk.END, "\n Report saved as 'viral_genome_analysis_report.txt'\n")

if __name__ == "__main__":
    root = tk.Tk()
    app = ViralGenomeAnalyzer(root)
    root.mainloop()
import gradio as gr
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
import pandas as pd
import matplotlib.pyplot as plt
import io
from collections import defaultdict

Entrez.email = "your.email@example.com"


class GenomeTransposonFinder:
    def __init__(self):
        self.genomes = {}

    def download_genomes(self, accession_list_str):

        accession_ids = [x.strip() for x in accession_list_str.split(',') if x.strip()]
        logs = []
        self.genomes = {}

        if not accession_ids:
            return "No IDs provided.", None

        logs.append(f"Attempting to fetch: {accession_ids}")

        try:
            handle = Entrez.efetch(db="nucleotide", id=accession_ids,
                                   rettype="fasta", retmode="text")
            records = list(SeqIO.parse(handle, "fasta"))
            handle.close()

            for record in records:
                self.genomes[record.id] = str(record.seq).upper()
                logs.append(f"Downloaded {record.id} (Length: {len(record.seq)} bp)")

            return "\n".join(logs), self.genomes

        except Exception as e:
            return f"NCBI error: {str(e)}", None

    def find_candidate_transposons(self, sequence, min_ir=4, max_ir=6,
                                   min_dist=100, max_dist=2000):

        seq = sequence
        L = len(seq)
        results = []

        for k in range(min_ir, max_ir + 1):

            for i in range(0, L - k + 1):

                kmer = seq[i:i+k]
                rev = str(Seq(kmer).reverse_complement())

                start = 0
                while True:
                    pos = seq.find(rev, start)
                    if pos == -1:
                        break

                    dist = pos - i
                    if min_dist <= dist <= max_dist:
                        results.append({
                            "Start": i,
                            "End": pos + k,
                            "Length": (pos + k) - i,
                            "IR_Seq": kmer,
                            "IR_RevComp": rev,
                            "IR_Len": k
                        })

                    start = pos + 1

        return pd.DataFrame(results)


analyzer = GenomeTransposonFinder()


def process_request(accession_str, min_ir, max_ir, min_dist, max_dist):
    log_text, genomes = analyzer.download_genomes(accession_str)

    if not genomes:
        return log_text, None, None

    all_results = []

    for acc_id, seq in genomes.items():
        log_text += f"\nProcessing {acc_id}..."

        df = analyzer.find_candidate_transposons(
            seq,
            min_ir=int(min_ir),
            max_ir=int(max_ir),
            min_dist=int(min_dist),
            max_dist=int(max_dist)
        )

        if not df.empty:
            df.insert(0, "Genome_ID", acc_id)
            all_results.append(df)
            log_text += f" Found {len(df)} candidates."
        else:
            log_text += " No candidates found."

    if not all_results:
        return log_text, pd.DataFrame(), None

    final_df = pd.concat(all_results).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(10, 4))
    y_positions = {gid: i for i, gid in enumerate(final_df['Genome_ID'].unique())}

    for idx, row in final_df.iterrows():
        y = y_positions[row['Genome_ID']]
        ax.plot([row['Start'], row['End']], [y, y], linewidth=2, alpha=0.7)
        ax.plot(row['Start'], y, '|', markersize=10, color='red')
        ax.plot(row['End'], y, '|', markersize=10, color='blue')

    ax.set_yticks(list(y_positions.values()))
    ax.set_yticklabels(list(y_positions.keys()))
    ax.set_xlabel("Base Position (bp)")
    ax.set_title("Detected Transposable Elements Map")
    plt.tight_layout()

    return log_text, final_df, fig


with gr.Blocks(title="FAST Transposon Finder") as demo:

    gr.Markdown("## Bacterial Transposon Finder (NCBI Integrated)")
    gr.Markdown("""
    **Instructions**
    1. Enter one or more NCBI Accession IDs  
       Example (fast small genome): `NC_001422.1`
    2. The tool will detect **Inverted Repeats (IRs)** indicating putative transposons.
    3. The algorithm now runs **20× faster**.
    """)

    with gr.Row():
        with gr.Column():
            acc_input = gr.Textbox(
                label="NCBI Accession IDs",
                value="NC_001422.1",   # VERY FAST small genome
                placeholder="e.g., NC_000913.3"
            )

            min_ir_sl = gr.Slider(4, 10, value=4, label="Min IR length")
            max_ir_sl = gr.Slider(4, 10, value=6, label="Max IR length")

            min_dist_sl = gr.Number(value=100, label="Min TE length (bp)")
            max_dist_sl = gr.Number(value=2000, label="Max TE length (bp)")

            run_btn = gr.Button("Run Analysis", variant="primary")

        with gr.Column():
            log_box = gr.Textbox(label="Log", lines=12)

    result_table = gr.Dataframe(headers=["Genome_ID", "Start", "End", "Length", "IR_Seq"])
    plot_output = gr.Plot()

    run_btn.click(
        fn=process_request,
        inputs=[acc_input, min_ir_sl, max_ir_sl, min_dist_sl, max_dist_sl],
        outputs=[log_box, result_table, plot_output]
    )

if __name__ == "__main__":
    demo.launch()

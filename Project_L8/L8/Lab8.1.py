import random
import gradio as gr


class TransposonSimulator:
    def __init__(self):
        self.bases = ['A', 'C', 'G', 'T']
        self.ir_seq = "GTCCGA"
        self.ir_seq_complement = "TCGGAC"
        self.tsd_length = 4

    def generate_random_dna(self, length):
        return "".join(random.choice(self.bases) for _ in range(length))

    def create_transposon_payload(self, length=20):
        return "".join(random.choice(self.bases) for _ in range(length))

    def insert_transposon(self, host_dna, index, label):
        if index + self.tsd_length > len(host_dna):
            return host_dna, f"Error: Index {index} out of bounds."

        target_site = host_dna[index: index + self.tsd_length]
        payload = self.create_transposon_payload()

        transposon_block = (
                self.ir_seq +
                payload +
                self.ir_seq_complement
        )

        new_dna = (
                host_dna[:index + self.tsd_length] +
                transposon_block +
                target_site +
                host_dna[index + self.tsd_length:]
        )

        log = (f"Inserted {label} at index {index}.\n"
               f"   Target Site (TSD): {target_site}\n"
               f"   Payload: {payload}")
        return new_dna, log

    def detect_transposons(self, dna_sequence):
        detected = []
        starts = []
        ends = []

        for i in range(len(dna_sequence)):
            if dna_sequence.startswith(self.ir_seq, i):
                starts.append(i)
            if dna_sequence.startswith(self.ir_seq_complement, i):
                ends.append(i + len(self.ir_seq_complement))

        for s in starts:
            for e in ends:
                if e > s:
                    tsd_start_idx = s - self.tsd_length
                    if tsd_start_idx < 0 or e + self.tsd_length > len(dna_sequence):
                        continue

                    tsd_left = dna_sequence[tsd_start_idx: s]
                    tsd_right = dna_sequence[e: e + self.tsd_length]

                    if tsd_left == tsd_right:
                        detected.append({
                            'start': s,
                            'end': e,
                            'tsd': tsd_left,
                            'len': e - s
                        })

        detected.sort(key=lambda x: x['len'])
        return detected


def run_simulation(dna_length, num_elements, do_intersection):
    sim = TransposonSimulator()
    logs = []

    dna = sim.generate_random_dna(int(dna_length))
    logs.append(f"Generated Host DNA: {len(dna)} bases.")

    current_dna = dna

    for i in range(int(num_elements)):
        current_len = len(current_dna)

        insert_idx = random.randint(10, int(current_len * 0.8))

        current_dna, log_msg = sim.insert_transposon(current_dna, insert_idx, f"Element_{i + 1}")
        logs.append(log_msg)

    if do_intersection:
        mid_point = len(current_dna) // 2
        current_dna, log_msg = sim.insert_transposon(current_dna, mid_point, "NESTED_ELEMENT")
        logs.append(f"FORCING INTERSECTION\n{log_msg}")

    detected = sim.detect_transposons(current_dna)

    result_text = f"Final DNA Sequence (Length {len(current_dna)}):\n{current_dna}\n\n"
    result_text += "=" * 30 + "\nDETECTION RESULTS\n" + "=" * 30 + "\n"

    if not detected:
        result_text += "No transposons detected."

    for i, item in enumerate(detected):
        seq_snippet = current_dna[item['start']:item['end']]
        display_seq = seq_snippet[:10] + "..." + seq_snippet[-10:] if len(seq_snippet) > 20 else seq_snippet

        result_text += (f"Found Transposon #{i + 1}\n"
                        f"  Position: {item['start']} - {item['end']}\n"
                        f"  Flanking TSD: {item['tsd']}\n"
                        f"  Content: {display_seq}\n"
                        f"  ----------------------------\n")

    return "\n".join(logs), result_text


with gr.Blocks(title="Transposon Simulator") as demo:
    gr.Markdown("# Transposon (Jumping Gene) Simulator")
    gr.Markdown(
        "Generates random DNA, inserts transposable elements with Inverted Repeats (IR) and Target Site Duplications (TSD), and attempts to detect them.")

    with gr.Row():
        with gr.Column():
            length_slider = gr.Slider(minimum=200, maximum=1000, value=300, step=50, label="Initial DNA Length")
            num_slider = gr.Slider(minimum=1, maximum=5, value=3, step=1, label="Number of Transposons")
            intersect_check = gr.Checkbox(label="Intersect (Nest) an Element?", value=True,
                                          info="Inserts an extra element inside another to test the algorithm (as per whiteboard note).")
            run_btn = gr.Button("Run Simulation", variant="primary")

        with gr.Column():
            log_output = gr.Textbox(label="Insertion Log", lines=10)

    result_output = gr.Textbox(label="Detection Results & DNA Sequence", lines=15, show_copy_button=True)

    run_btn.click(fn=run_simulation, inputs=[length_slider, num_slider, intersect_check],
                  outputs=[log_output, result_output])

if __name__ == "__main__":
    demo.launch()
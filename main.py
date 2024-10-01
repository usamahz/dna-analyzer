# DNA Sequence Analysis
# Author: Usamah Zaheer

from dna_analysis import load_sequences, analyse_sequences

if __name__ == "__main__":
    sequences = load_sequences("dna_sequences.json")
    summary = analyse_sequences(sequences)
    print(summary)

    # Save the summary
    with open("analysis_summary.md", "w") as f:
        f.write(summary)

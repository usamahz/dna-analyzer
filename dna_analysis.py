# DNA Sequence Analysis
# Author: Usamah Zaheer

import json
from collections import Counter
import logging

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def load_sequences(file_path):
    """
    Load DNA sequences from a JSON file.

    Args:
        file_path (str): Path to the JSON file containing DNA sequences.

    Returns:
        list: A list of DNA sequences.
    """
    with open(file_path, "r") as f:
        data = json.load(f)
    return data["sequences"]


def calculate_gc_content(sequence):
    """
    Calculate the GC content of a DNA sequence.

    Args:
        sequence (str): A DNA sequence.

    Returns:
        float: The GC content as a fraction of the total sequence length.
    """
    gc_count = sequence.count("G") + sequence.count("C")
    return gc_count / len(sequence)


def calculate_dinucleotide_frequencies(sequence):
    """
    Calculate the frequencies of dinucleotides in a DNA sequence.

    Args:
        sequence (str): A DNA sequence.

    Returns:
        Counter: A Counter object containing the frequencies of dinucleotides.
    """
    dinucleotides = [sequence[i : i + 2] for i in range(len(sequence) - 1)]
    return Counter(dinucleotides)


def find_most_common_kmers(sequence, k):
    """
    Find the most common k-mers in a DNA sequence.

    Args:
        sequence (str): A DNA sequence.
        k (int): The length of the k-mers to find.

    Returns:
        list: A list of tuples containing the 5 most common k-mers and their frequencies.
    """
    kmers = [sequence[i : i + k] for i in range(len(sequence) - k + 1)]
    return Counter(kmers).most_common(5)


def find_palindromes(sequence, min_length=20):
    palindromes = []
    for i in range(len(sequence)):
        # Odd-length palindromes
        left, right = i, i
        while left >= 0 and right < len(sequence) and sequence[left] == sequence[right]:
            if (
                right - left + 1 >= min_length
                and len(set(sequence[left : right + 1])) >= 3
            ):
                palindromes.append(sequence[left : right + 1])
            left -= 1
            right += 1

        # Even-length palindromes
        left, right = i, i + 1
        while left >= 0 and right < len(sequence) and sequence[left] == sequence[right]:
            if (
                right - left + 1 >= min_length
                and len(set(sequence[left : right + 1])) >= 3
            ):
                palindromes.append(sequence[left : right + 1])
            left -= 1
            right += 1

    return palindromes


def analyse_sequences(sequences):
    """
    Analyse a list of DNA sequences and generate a summary of various statistics and patterns.

    Args:
        sequences (list): A list of DNA sequences to analyse.

    Returns:
        str: A formatted string containing the analysis summary.
    """
    gc_contents = []
    overall_dinucleotide_freq = Counter()
    kmers = {k: Counter() for k in [3, 4, 5]}
    palindromes = []
    sequence_lengths = []
    nucleotide_composition = Counter()

    logging.info("Analyzing sequences...")
    for i, seq in enumerate(sequences):
        logging.debug(f"Processing sequence {i+1}/{len(sequences)}")
        # GC content
        gc_count = seq.count("G") + seq.count("C")
        gc_contents.append(gc_count / len(seq))

        # Dinucleotide frequencies
        overall_dinucleotide_freq.update(seq[i : i + 2] for i in range(len(seq) - 1))

        # k-mers
        for k in kmers:
            kmers[k].update(seq[i : i + k] for i in range(len(seq) - k + 1))

        # Palindromes
        palindromes.extend((i, p) for p in find_palindromes(seq))

        # Sequence length
        sequence_lengths.append(len(seq))

        # Nucleotide composition
        nucleotide_composition.update(seq)

    avg_length = sum(sequence_lengths) / len(sequences)
    min_length = min(sequence_lengths)
    max_length = max(sequence_lengths)

    # Modify the sequence length part of the summary
    if min_length == max_length:
        length_summary = f"""
       a. Sequence Length:
          All sequences have the same length: {min_length} BP"""
    else:
        length_summary = f"""
       a. Sequence Length:
          - Average length: {avg_length:.2f} BP
          - Shortest sequence: {min_length} BP
          - Longest sequence: {max_length} BP"""

    repeat_sequences = find_repeat_sequences(sequences)

    logging.info("Generating summary.....")
    # Generate summary
    summary = f"""
    Learning Model of Life: 200 DNA sequences Analysis Summary

    1. GC Content Distribution:
       - Mean: {sum(gc_contents) / len(gc_contents):.2f}
       - Min: {min(gc_contents):.2f}
       - Max: {max(gc_contents):.2f}
       - Distribution:
         {' '.join(f'{gc:.2f}: {count:3d}' for gc, count in sorted(Counter(min([round(i * 0.05, 2) for i in range(21)], key=lambda x: abs(x - gc)) for gc in gc_contents).items()))}
    
    2. Dinucleotide Frequencies:
       {dict(overall_dinucleotide_freq)}
    
    3. 5 Most Common k-mers:
       3-mers: {kmers[3].most_common(5)}
       4-mers: {kmers[4].most_common(5)}
       5-mers: {kmers[5].most_common(5)}
    
    4. Palindromes, above 20 base pairs in length that include at least 3 of the four bases (A,C,T,G):
       Total found: {len(palindromes)}
       All palindromes: {palindromes}
    
    5. Additional Observations:
       {length_summary}
       
       b. Repeat Sequences:
          - Total unique repeat sequences found (of length 20): {len(repeat_sequences)}
          - Top 10:
            {repeat_sequences.most_common(10)}
       
       c. Nucleotide Composition:
          {dict(nucleotide_composition)}
    """

    logging.info("Analysis complete.")
    return summary


def find_repeat_sequences(sequences, min_length=20):
    """Find the top 5 most frequent repeat sequences across all input sequences."""
    repeat_counts = Counter()

    logging.info("Finding repeat sequences...")
    for seq in sequences:
        for i in range(len(seq) - min_length + 1):
            substr = seq[i : i + min_length]
            repeat_counts[substr] += 1

    return repeat_counts


def calculate_nucleotide_composition(sequences):
    """Calculate the overall nucleotide composition of all sequences."""
    total_count = sum(len(seq) for seq in sequences)
    composition = Counter("".join(sequences))
    return {base: count / total_count for base, count in composition.items()}
# DNA Sequence Analysis

## Author: Usamah Zaheer

This project performs an analysis of given DNA sequences, providing various statistics and patterns found within the genetic data.

## Project Structure

- `main.py`: The entry point of the program.
- `dna_analysis.py`: Contains the core functions for DNA sequence analysis.
- `dna_sequences.json`: Input file containing the DNA sequences to be analysed.
- `analysis_summary.md`: Output file containing the analysis results.

## Features

The analysis is segregated into the following sections:

a. Basic sequence statistics:
   - Overall GC content distribution
   - Dinucleotide frequencies

b. Identification of the top 5 most common k-mers (substrings) for k=3, 4, and 5

c. Detection of unusual patterns, such as:
   - Palindromic sequences above 20 base pairs in length that include at least 3 of the four bases (A,C,T,G)
   - Any other patterns found interesting or unusual

d. A brief summary of findings, highlighting any sequences or patterns that stand out.

## How to Run

1. Ensure you have Python 3.x installed on your system.
2. Place your DNA sequences in the `dna_sequences.json` file.
3. Run the following command in the terminal:

   ```
   python main.py
   ```

4. The analysis results will be printed to the console and saved in `analysis_summary.md`.

## Dependencies

This project uses only Python standard libraries and does not require any additional installations.

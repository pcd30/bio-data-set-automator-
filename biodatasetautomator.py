# -*- coding: utf-8 -*-
"""BioDatasetAutomator



Original file is located at
    https://colab.research.google.com/drive/14O_Zc1t2HubqUrrag7wfK0pK1B42I55a
"""

!pip install biopython pandas matplotlib

import random

nucleotides = ["A","T","G","C"]

def generate_sequence(length):
    return "".join(random.choice(nucleotides) for _ in range(length))

for i in range(5):
    seq = generate_sequence(50)
    print(f">Sequence_{i+1}")
    print(seq)

import random

nucleotides = ["A","T","G","C"]

def generate_sequence(length):
    return "".join(random.choice(nucleotides) for _ in range(length))

with open("random_sequences.fasta", "w") as file:
    for i in range(20):
        seq = generate_sequence(100)
        file.write(f">Sequence_{i+1}\n")
        file.write(seq + "\n")

print("FASTA file created")

from Bio import SeqIO

for record in SeqIO.parse("random_sequences.fasta", "fasta"):
    print(record.id, len(record.seq))

def gc_content(seq):
    g = seq.count("G")
    c = seq.count("C")
    return (g+c)/len(seq)*100

from Bio import SeqIO

for record in SeqIO.parse("random_sequences.fasta", "fasta"):
    gc = gc_content(record.seq)
    print(record.id, "GC:", round(gc,2),"%")

import matplotlib.pyplot as plt
from Bio import SeqIO

gc_values = []

def gc_content(seq):
    g = seq.count("G")
    c = seq.count("C")
    return (g+c)/len(seq)*100

for record in SeqIO.parse("random_sequences.fasta", "fasta"):
    gc_values.append(gc_content(record.seq))

plt.hist(gc_values)
plt.xlabel("GC Content")
plt.ylabel("Number of Sequences")
plt.title("GC Content Distribution")
plt.show()

from Bio import SeqIO

motif = "ATG"

for record in SeqIO.parse("random_sequences.fasta", "fasta"):
    sequence = str(record.seq)
    positions = []

    for i in range(len(sequence) - len(motif) + 1):
        if sequence[i:i+len(motif)] == motif:
            positions.append(i)

    print(record.id, "Motif positions:", positions)

from Bio import SeqIO

motif = "ATG"

for record in SeqIO.parse("random_sequences.fasta", "fasta"):
    sequence = str(record.seq)
    count = sequence.count(motif)

    print(record.id, "ATG count:", count)

from Bio import SeqIO
import numpy as np

lengths = []

for record in SeqIO.parse("random_sequences.fasta", "fasta"):
    lengths.append(len(record.seq))

print("Total sequences:", len(lengths))
print("Average length:", np.mean(lengths))
print("Maximum length:", max(lengths))
print("Minimum length:", min(lengths))

import pandas as pd
from Bio import SeqIO

def gc_content(seq):
    g = seq.count("G")
    c = seq.count("C")
    return (g+c)/len(seq)*100

data = []

for record in SeqIO.parse("random_sequences.fasta", "fasta"):
    data.append({
        "Sequence_ID": record.id,
        "Length": len(record.seq),
        "GC_Content": gc_content(record.seq)
    })

df = pd.DataFrame(data)

df

df.to_csv("sequence_analysis_results.csv", index=False)

print("Results saved.")

from google.colab import files
files.download("sequence_analysis_results.csv")

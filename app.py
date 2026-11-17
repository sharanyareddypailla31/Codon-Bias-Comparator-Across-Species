import io
import math
from collections import Counter, defaultdict
from itertools import combinations

import pandas as pd
import plotly.express as px
import streamlit as st
from Bio import SeqIO

st.set_page_config(page_title="Codon Bias Comparator", page_icon="🧬", layout="wide")
st.title("🧬 Codon Bias Comparator Across Species")
st.write("Compare codon usage in FASTA files from 2–4 organisms.")
st.info("Upload coding DNA sequences (CDS). Each FASTA record is analysed separately.")

# Standard genetic code
CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L", "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*", "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L", "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q", "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M", "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K", "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V", "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E", "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

# Group codons that code for the same amino acid (stop codons excluded).
SYNONYMOUS = defaultdict(list)
for codon, amino_acid in CODON_TABLE.items():
    if amino_acid != "*":
        SYNONYMOUS[amino_acid].append(codon)


def read_fasta(file):
    """Read and check the uploaded FASTA file."""
    records = list(SeqIO.parse(io.StringIO(file.getvalue().decode("utf-8")), "fasta"))
    if not records:
        raise ValueError("No FASTA records found. Check the file format.")

    sequences = []
    for record in records:
        sequence = str(record.seq).upper().replace("U", "T")
        if set(sequence) - set("ACGT"):
            raise ValueError(f"Invalid DNA letters in record: {record.id}")
        if len(sequence) < 3:
            raise ValueError(f"Record {record.id} is shorter than one codon.")
        sequences.append(sequence)
    return sequences

  
def analyse(sequences):
    """Calculate codon counts, frequencies, GC content, GC3 and RSCU."""
    counts = Counter()
    total_bases = gc_bases = gc3_bases = codon_total = trailing = 0

    for sequence in sequences:
        total_bases += len(sequence)
        gc_bases += sequence.count("G") + sequence.count("C")
        usable_length = len(sequence) - len(sequence) % 3
        trailing += len(sequence) % 3

        for i in range(0, usable_length, 3):
            codon = sequence[i:i + 3]
            counts[codon] += 1
            codon_total += 1
            if codon[2] in "GC":
                gc3_bases += 1

    frequencies = {
        codon: counts[codon] / codon_total * 100
        for codon in CODON_TABLE
    }

    rscu = {}
    for codons in SYNONYMOUS.values():
        family_total = sum(counts[c] for c in codons)
        for codon in codons:
            rscu[codon] = counts[codon] * len(codons) / family_total if family_total else 0

    return {
        "counts": counts,
        "frequencies": frequencies,
        "rscu": rscu,
        "gc": gc_bases / total_bases * 100,
        "gc3": gc3_bases / codon_total * 100 if codon_total else 0,
        "length": total_bases,
        "codons": codon_total,
        "trailing": trailing,
        "records": len(sequences),
    }


def similarity(a, b):
    """Custom cosine similarity of synonymous-codon count patterns."""
    scores = []
    for codons in SYNONYMOUS.values():
        x = [a["counts"][c] for c in codons]
        y = [b["counts"][c] for c in codons]
        nx = math.sqrt(sum(v * v for v in x))
        ny = math.sqrt(sum(v * v for v in y))
        if nx and ny:
            scores.append(sum(i * j for i, j in zip(x, y)) / (nx * ny))
    return sum(scores) / len(scores) * 100 if scores else 0


# Choose organisms and collect their FASTA files.
number = st.selectbox("Number of organisms", [2, 3, 4])
default_names = ["Human", "E. coli", "Yeast", "Arabidopsis"]
names, uploads = [], []
columns = st.columns(2)

for i in range(number):
    with columns[i % 2]:
        name = st.text_input(f"Organism {i + 1}", default_names[i], key=f"name{i}").strip()
        upload = st.file_uploader(f"FASTA file for {name or 'organism'}", ["fasta", "fa", "fna", "txt"], key=f"file{i}")
        names.append(name or f"Organism {i + 1}")
        uploads.append(upload)

if not all(uploads):
    st.warning("Upload one FASTA file for every selected organism.")
elif len(set(names)) != len(names):
    st.error("Please use a different name for each organism.")
else:
    try:
        data = {
            name: analyse(read_fasta(file))
            for name, file in zip(names, uploads)
        }
        st.success("All FASTA files analysed successfully.")

        # 1. Sequence summary
        st.header("1. Sequence summary")
        summary = pd.DataFrame([
            {
                "Organism": name,
                "FASTA records": result["records"],
                "Sequence length (nt)": result["length"],
                "Complete codons": result["codons"],
                "GC content (%)": round(result["gc"], 2),
                "GC3 (%)": round(result["gc3"], 2),
                "Unused trailing bases": result["trailing"],
            }
            for name, result in data.items()
        ])
        st.dataframe(summary, use_container_width=True, hide_index=True)

        # 2. Codon usage
        st.header("2. Codon usage comparison")
        codon_rows = []
        for codon, aa in CODON_TABLE.items():
            row = {"Codon": codon, "Amino acid": "Stop" if aa == "*" else aa}
            for name in names:
                row[f"{name} count"] = data[name]["counts"][codon]
                row[f"{name} frequency (%)"] = round(data[name]["frequencies"][codon], 3)
            codon_rows.append(row)

        codon_df = pd.DataFrame(codon_rows)
        st.dataframe(codon_df, use_container_width=True, hide_index=True)
        st.download_button("Download codon comparison CSV", codon_df.to_csv(index=False), "codon_comparison.csv")

        frequency_df = pd.DataFrame([
            {"Codon": codon, "Organism": name, "Frequency (%)": data[name]["frequencies"][codon]}
            for name in names for codon in CODON_TABLE
        ])
        st.plotly_chart(px.bar(
            frequency_df, x="Codon", y="Frequency (%)", color="Organism",
            barmode="group", title="Codon frequency by organism"
        ), use_container_width=True)

        # 3. RSCU
        st.header("3. Relative Synonymous Codon Usage (RSCU)")
        rscu_df = pd.DataFrame([
            {"Codon": codon, "Amino acid": aa,
             **{name: data[name]["rscu"][codon] for name in names}}
            for aa, codons in SYNONYMOUS.items() for codon in codons
        ])
        st.dataframe(rscu_df.round(3), use_container_width=True, hide_index=True)
        st.caption("RSCU > 1: used more than expected; RSCU < 1: used less than expected; 0: not observed.")
        st.download_button("Download RSCU CSV", rscu_df.to_csv(index=False), "rscu_comparison.csv")
        st.plotly_chart(px.imshow(
            rscu_df.set_index("Codon")[names], aspect="auto",
            color_continuous_scale="RdBu_r", title="RSCU heatmap"
        ), use_container_width=True)

        # 4. Codons with the biggest frequency differences
        st.header("4. Most different codon preferences")
        differences = []
        for codon, aa in CODON_TABLE.items():
            if aa == "*":
                continue
            values = {name: data[name]["frequencies"][codon] for name in names}
            high, low = max(values, key=values.get), min(values, key=values.get)
            differences.append({
                "Codon": codon, "Amino acid": aa,
                "Highest frequency organism": high,
                "Highest frequency (%)": round(values[high], 3),
                "Lowest frequency organism": low,
                "Lowest frequency (%)": round(values[low], 3),
                "Frequency range (%)": round(values[high] - values[low], 3),
            })
        differences_df = pd.DataFrame(differences).sort_values("Frequency range (%)", ascending=False)
        st.dataframe(differences_df.head(10), use_container_width=True, hide_index=True)

        # 5. Pairwise similarity and interpretation
        st.header("5. Overall similarity and interpretation")
        pairs = [(a, b, similarity(data[a], data[b])) for a, b in combinations(names, 2)]
        matrix = pd.DataFrame(100.0, index=names, columns=names)
        for a, b, score in pairs:
            matrix.loc[a, b] = matrix.loc[b, a] = score

        st.subheader("Codon preference similarity (%)")
        st.dataframe(matrix.round(2), use_container_width=True)
        st.plotly_chart(px.imshow(
            matrix, text_auto=".2f", zmin=0, zmax=100, aspect="auto",
            color_continuous_scale="RdYlGn", title="Pairwise similarity heatmap"
        ), use_container_width=True)

        highest = max(pairs, key=lambda item: item[2])
        lowest = min(pairs, key=lambda item: item[2])
        st.write(f"**Most similar:** {highest[0]} and {highest[1]} ({highest[2]:.2f}%).")
        st.write(f"**Least similar:** {lowest[0]} and {lowest[1]} ({lowest[2]:.2f}%).")

        st.markdown("**Interpretation**")
        for a, b, score in sorted(pairs, key=lambda item: item[2], reverse=True):
            meaning = ("high" if score >= 80 else
                       "moderate-to-high" if score >= 60 else
                       "moderate" if score >= 40 else "relatively low")
            st.write(f"- {a} and {b}: {score:.2f}% ({meaning} similarity).")

        st.warning(
            "This is a custom cosine-based codon preference score, not CAI or DNA identity. "
            "Results depend on the uploaded genes and sequence quality."
        )

    except Exception as error:
        st.error(f"Could not analyse the files: {error}")

st.divider()
st.caption("Use comparable coding sequences in the correct reading frame for meaningful biological comparisons.")

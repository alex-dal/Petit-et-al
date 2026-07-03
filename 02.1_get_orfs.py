#!/usr/bin/env python3
from Bio import SeqIO

# Parameters
#input_fasta = "PcamVGP1_loci_flanks.fasta"
#output_fasta = "PcamVGP1_loci_flanks_300longorfs.faa"
input_fasta = "PcamVGP1_starships.fasta"
output_fasta = "PcamVGP1_starships_100longorfs.faa"

min_aa = 100  # minimum ORF length in amino acids

with open(output_fasta, "w") as out_f:
    for record in SeqIO.parse(input_fasta, "fasta"):
        seq = record.seq
        for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
            for frame in range(3):
                seq_trimmed = nuc[frame:]
                # trim last 1-2 bases to make length divisible by 3
                seq_trimmed = seq_trimmed[:len(seq_trimmed) - (len(seq_trimmed) % 3)]
                trans = seq_trimmed.translate(to_stop=False)
                
                aa_start = 0
                for i, aa in enumerate(trans):
                    if aa == "*":
                        if i - aa_start >= min_aa:
                            # Calculate DNA coordinates
                            start = frame + aa_start*3 + 1
                            end = frame + i*3
                            if strand == -1:
                                start, end = len(seq) - end + 1, len(seq) - start + 1
                            # Write FASTA
                            out_f.write(f">{record.id}_strand{strand}_frame{frame+1}_{start}_{end}\n")
                            out_f.write(str(trans[aa_start:i]) + "\n")
                        aa_start = i + 1
                # handle ORF at end without stop codon
                if len(trans) - aa_start >= min_aa:
                    start = frame + aa_start*3 + 1
                    end = frame + len(trans)*3
                    if strand == -1:
                        start, end = len(seq) - end + 1, len(seq) - start + 1
                    out_f.write(f">{record.id}_strand{strand}_frame{frame+1}_{start}_{end}\n")
                    out_f.write(str(trans[aa_start:]) + "\n")

##############################################################################
# Extract sequences of Starships that contain PcamV-GP1 insertions
##############################################################################

# Extract sequences
bedtools getfasta -fi all_genomes_combined.fna -bed PcamVGP1_starships.bed -fo PcamVGP1_starships.fasta -name -s

# Translate filtered blast hits into 6 frames, keep 300aa+
python3 get_orfs.py

# Split this output in batches of 999
seqkit split -i -s 999 PcamVGP1_starships_100longorfs.faa -O cdsearch_batches_starships/
# Remove #-ed header lines

# Clean cdsearch outputs
for file in cdsearch_batches_starships/*.txt; do
python3 clean_cdsearch_outputs.py $file ;
done

# Merge cdsearch outputs, keep one header
(first_file=$(ls cdsearch_batches_starships/*_cln.txt | head -n1); head -n1 "$first_file"; for f in cdsearch_batches_starships/*_cln.txt; do tail -n +2 "$f"; done) > starship_cdsearch_outs_cln.txt

# Merge identical protein domains if they overlap, keep longest domain if there's an overlap, and compute DNA coordinates of the protein domains
python3 merge_cdsearch_domains.py

# Format as bed [s=bed_Start e=bed_end]
awk -F '\t' 'BEGIN{OFS="\t"}
NR>1 && $14!="" && $15!="" {

    strand = ($12==1 ? "+" : "-")

    s = ($16 < $17 ? $16 : $17)
    e = ($16 > $17 ? $16 : $17)

    print $1, s-1, e, $9, $7, strand
}' starship_cdsearch_outs_merged_filtered.txt > starship_cdsearch_outs_merged_filtered.bed


# Simplify headers of old fasta file
sed '/^>/ s/^>.*::/>/' PcamVGP1_starships.fasta > PcamVGP1_starships_fixhead.fasta
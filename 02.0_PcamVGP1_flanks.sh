
##############################################################################
# Extract 50kb upstream and downstream of viral loci
##############################################################################

samtools faidx all_genomes_combined.fna
cut -f1,2 all_genomes_combined.fna.fai > all_genomes_combined.genome

# Make flank bed
bedtools slop \
    -i hits.bed \
    -g all_genomes_combined.fna.fai \
    -b 50000 \
> PcamVGP1_loci_flanks.bed

# Extract sequences
bedtools getfasta -fi all_genomes_combined.fna -bed PcamVGP1_loci_flanks.bed -fo PcamVGP1_loci_flanks.fasta -name -s

# Translate filtered blast hits into 6 frames, keep 300aa+
python3 get_orfs.py

# Split this output in batches of 999
seqkit split -i -s 999 PcamVGP1_loci_flanks_300longorfs.faa -O cdsearch_batches/

# Clean cdsearch outputs
for file in cdsearch_batches/*.txt; do
python3 clean_cdsearch_outputs.py $file ;
done

# Merge cdsearch outputs, keep one header
(first_file=$(ls cdsearch_batches/*_cln.txt | head -n1); head -n1 "$first_file"; for f in cdsearch_batches/*_cln.txt; do tail -n +2 "$f"; done) > cdsearch_outs_cln.txt

# Merge identical protein domains if they overlap, keep longest domain if there's an overlap, and compute DNA coordinates of the protein domains
python3 merge_cdsearch_domains.py

# Format as bed [s=bed_Start e=bed_end]
awk -F '\t' 'BEGIN{OFS="\t"}
NR>1 && $14!="" && $15!="" {

    strand = ($12==1 ? "+" : "-")

    s = ($16 < $17 ? $16 : $17)
    e = ($16 > $17 ? $16 : $17)

    print $1, s-1, e, $9, $7, strand
}' cdsearch_outs_merged.txt > cdsearch_outs_merged.bed


# Simplify headers of old fasta file
sed '/^>/ s/^>.*::/>/' PcamVGP1_loci_flanks.fasta > PcamVGP1_loci_flanks_fixhead.fasta
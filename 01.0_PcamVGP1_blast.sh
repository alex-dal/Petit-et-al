##############################################################################
# Prep files for blast
############################################################################### 

# Add prefix to non ncbi files
for dir in [0-9]*; do
    # Skip GCA/GCF directories
    if [[ "$dir" == GCA* ]] || [[ "$dir" == GCF* ]]; then
        continue
    fi

    # Loop over .fna files in the directory
    for f in "$dir"/*.fna; do
        # Skip if no .fna file exists
        [[ -e "$f" ]] || continue

        # Output file: add _prefixed to original filename
        out="${f%.fna}_prefixed.fna"

        # Add directory name as prefix to each header
        awk -v prefix="$dir" '/^>/{print ">"prefix"_"substr($0,2); next}{print}' "$f" > "$out"

        echo "Processed $f -> $out"
    done
done

# cat all genomes together
cat */*.fna > all_genomes_combined.fna

# Replace spaces in headers by _ and clean : , ( )
awk '
{
    if(/^>/){
        gsub(/[,:()]/,"")   # remove , : ( )
        gsub(/ /,"_")       # replace spaces with underscores
    }
    print
}' all_genomes_combined.fna > all_genomes_combined_cln.fna

# Check if headers are unique (some non-NCBI assemblies here)
awk '/^>/ {
    if($0 in seen) {print "Duplicate header: "$0}
    else {seen[$0]=1}
}'  all_genomes_combined_cln.fna

# They are.
mv all_genomes_combined_cln.fna all_genomes_combined.fna


############################################################################### 
# Blastn against genome db
############################################################################### 

# Blast db
makeblastdb -in all_genomes_combined.fna -dbtype nucl -out all_genomes_db

# Align
blastn \
  -query PcamV-GP1.fa \
  -db all_genomes_db \
  -evalue 1e-5 \
  -num_threads 8 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
  -out PcamV-GP1_vs_all_genomes.txt

# Keep only alignments >2kb
awk '$4 > 2000' PcamV-GP1_vs_all_genomes.txt > PcamV-GP1_vs_all_genomes_2000.txt

# Make bed file with ranges of hit loci
awk '{
    s=$9; e=$10;
    if(s>e){tmp=s;s=e;e=tmp; strand="-"} else strand="+";
    print $2"\t"(s-1)"\t"e"\t"$2"\t0\t"strand
}' PcamV-GP1_vs_all_genomes_2000.txt > hits.bed

# Bed file manually edited for early-diverging element: hits.bed2
# I extended the genomic regions of blast hits (start+6434 bp, respective of orientation), and this recovered protein domains that were cut off in the initial blast hits

# Get sequences, taking into account reverse alignments (-s)
bedtools getfasta -fi all_genomes_combined.fna -bed hits.bed -fo PcamV-GP1_hits_extracted.fna -name -s

# Manually paste PcamV-GP1 sequence to the fasta file


##############################################################################
# Align
##############################################################################

mafft --adjustdirection  PcamV-GP1_hits_extracted.fna >   PcamV-GP1_hits_extracted.aln-mafft


##############################################################################
# Construct maximum-likelihood trees
##############################################################################

iqtree3 -s PcamV-GP1_hits_extracted.aln-mafft -bb 1000 -alrt 1000 -nt AUTO


##############################################################################
# Annotate protein domains of fasta sequences
##############################################################################

# Translate filtered blast hits into 6 frames
seqkit translate -w 0 -f 6 PcamV-GP1_hits_extracted.fna > PcamV-GP1_hits_extracted.faa

# Keep only long ORFs with > 100 aa. If a sequence has more than 1 ORF >100aa, all are kept
bioawk -c fastx '
{
    split($seq, orfs, "\\*");
    for(i in orfs){
        if(length(orfs[i])>=100){
            print ">"$name"_orf"i"\n"orfs[i]
        }
    }
}' PcamV-GP1_hits_extracted.faa > PcamV-GP1_hits_extracted_100longorfs.faa

# Batch search for protein domains on https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi? ; default mode

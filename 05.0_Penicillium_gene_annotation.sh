
PROTEINS="~/protein_evidence.faa"

file=$(sed -n "${SLURM_ARRAY_TASK_ID}p" genomes_augustus.txt)

dir=$(dirname "$file")
filename=$(basename "$file")
base=${filename%.rm.fa}

paf="${dir}/${base}.miniprot.paf"
gff="${dir}/${base}.miniprot.gff"
augustus_gff="${dir}/${base}.augustus.gff"

echo "Processing: $file"

# for paf
miniprot -t 2 "$file" "$PROTEINS" > "$paf"

# for hints
awk 'BEGIN{OFS="\t"}
{
    if ($0 ~ /^#/ || NF < 12) next;

    chr=$1;
    start=$8+1;
    end=$9;
    strand=($5=="+") ? "+" : "-";

    print chr, "P", "CDSpart", start, end, ".", strand, ".", "source=P";
}' "$paf" > "$gff"

augustus \
    --species=anidulans \
    --softmasking=1 \
    --hintsfile="$gff" \
    --extrinsicCfgFile=$AUGUSTUS_CONFIG_PATH/extrinsic/extrinsic.MPE.cfg \
    "$file" \
    > "$augustus_gff"

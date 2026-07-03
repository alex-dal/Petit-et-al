==========================================================================================
---------------------------------------WORKFLOW-------------------------------------------
==========================================================================================

--------------SETUP--------------

conda activate starfish

STARFISHDIR=$(dirname -- $(command -v -- starfish))/../

realpath assembly/* | perl -pe 's/^(.+?([^\/]+?).fna)$/\2\t\1/' > ome2assembly.txt
realpath gff/* | perl -ne 'chomp; /\/([^\/]+)\.gff$/; print "$1\t$_\n"' > ome2gff.txt 

cat gff/*.gff > penicillium.gff //// GGF > GFF3

mkdir blastdb
cut -f2 ome2assembly.txt | xargs cat > blastdb/penicillium.assemblies.fna
makeblastdb -in blastdb/penicillium.assemblies.fna -out blastdb/penicillium.assemblies -parse_seqids -dbtype nucl

$STARFISHDIR/aux/seq-gc.sh -Nbw 1000 blastdb/penicillium.assemblies.fna > penicillium.assemblies.gcContent_w1000.bed
rm blastdb/penicillium.assemblies.fna

--------------GENE FINDER--------------

mkdir geneFinder


//LONG ASS STEP THIS TOOK 1H20MINS
starfish annotate -T 2 -x penicillium_tyr -a ome2assembly.txt -g ome2gff.txt -p $STARFISHDIR/db/YRsuperfams.p1-512.hmm -P $STARFISHDIR/db/YRsuperfamRefs.faa -i tyr -o geneFinder/

starfish consolidate -o ./ -g penicillium.gff -G geneFinder/penicillium_tyr.filt_intersect.gff

realpath penicillium_tyr.filt_intersect.consolidated.gff | perl -pe 's/^/macpha6\t/' > ome2consolidatedGFF.txt

starfish sketch -m 10000 -q geneFinder/penicillium_tyr.filt_intersect.ids -g ome2consolidatedGFF.txt -i s -x macpha6 -o geneFinder/

grep -P '\ttyr\t' geneFinder/penicillium.bed > geneFinder/penicillium.tyr.bed 

--------------ELEMENT FINDER--------------

mkdir elementFinder

//LONG ASS STEP THIS TOOK 1H
starfish insert -T 2 -a ome2assembly.txt -d blastdb/penicillium.assemblies -b geneFinder/penicillium.tyr.bed -i tyr -x macpha6 -o elementFinder/

//LONG ASS STEP THIS TOOK 1H30
starfish insert -T 2 -a ome2assembly.txt -d blastdb/penicillium.assemblies -b elementFinder/penicillium.insert.bed -i tyr -x penicillium_round2 -o elementFinder/ --pid 80 --hsp 750

starfish summarize -a ome2assembly.txt -b elementFinder/penicillium.insert.bed -x macpha6 -o elementFinder/ -S elementFinder/penicillium.insert.stats -g ome2consolidatedGFF.txt -t geneFinder/penicillium_tyr.filt_intersect.ids

awk '{
    if ($4 ~ /DR/) print $0 "\t255,255,0";
    else if ($4 ~ /TIR/) print $0 "\t255,165,0";
    else if ($5 ~ /tyr|cap/) print $0 "\t255,0,0";
    else if ($5 != ".") print $0 "\t128,0,128";
    else print $0 "\t169,169,169";
}' elementFinder/penicillium.elements.bed > elementFinder/penicillium.elements.color.bed

hmmsearch --noali --notextw -E 0.001 --max --cpu 12 --tblout elementFinder/penicillium_tyr_vs_YRsuperfams.out $STARFISHDIR/db/YRsuperfams.p1-512.hmm geneFinder/penicillium_tyr.filt_intersect.fas
perl -p -e 's/ +/\t/g' elementFinder/penicillium_tyr_vs_YRsuperfams.out | cut -f1,3,5 | grep -v '#' | sort -k3,3g | awk '!x[$1]++' > elementFinder/penicillium_tyr_vs_YRsuperfams_besthits.txt

grep -P '\tcap\t' elementFinder/penicillium.elements.bed | cut -f4,7 > elementFinder/penicillium.cap2ship.txt
$STARFISHDIR/aux/searchReplace.pl --strict -i elementFinder/penicillium_tyr_vs_YRsuperfams_besthits.txt -r elementFinder/penicillium.cap2ship.txt > elementFinder/penicillium_elements_vs_YRsuperfams_besthits.txt

mmseqs easy-cluster geneFinder/penicillium_tyr.filt_intersect.fas elementFinder/penicillium_tyr elementFinder/ --threads 2 --min-seq-id 0.5 -c 0.25 --alignment-mode 3 --cov-mode 0 --cluster-reassign
$STARFISHDIR/aux/mmseqs2mclFormat.pl -i elementFinder/penicillium_tyr_cluster.tsv -g navis -o elementFinder/

starfish sim -m element -t nucl -b elementFinder/penicillium.elements.bed -x macpha6 -o elementFinder/ -a ome2assembly.txt
starfish group -m mcl -s elementFinder/penicillium.element.nucl.sim -i hap -o elementFinder/ -t 0.05

$STARFISHDIR/aux/searchReplace.pl -i elementFinder/penicillium_tyr_cluster.mcl -r elementFinder/penicillium.cap2ship.txt > elementFinder/penicillium.element_cluster.mcl

$STARFISHDIR/aux/mergeGroupfiles.pl -t elementFinder/penicillium.element_cluster.mcl -q elementFinder/penicillium.element.nucl.I1.5.mcl > elementFinder/penicillium.element.navis-hap.mcl

awk '{ for (i = 2; i <= NF; i++) print $i"\t"$1 }' elementFinder/penicillium.element.navis-hap.mcl > elementFinder/penicillium.element.navis-hap.txt

join -t$'\t' -1 1 -2 2 <(sort -t$'\t' -k1,1 elementFinder/penicillium.element.navis-hap.txt | grep -P '_e|_s') <(sort -t$'\t' -k2,2 elementFinder/penicillium.elements.feat) | awk -F'\t' '{print}' > elementFinder/penicillium.elements.temp.feat
echo -e "#elementID\tfamilyID\tnavisHapID\tcontigID\tcaptainID\telementBegin\telementEnd\telementLength\tstrand\tboundaryType\temptySiteID\temptyContig\temptyBegin\temptyEnd\temptySeq\tupDR\tdownDR\tDRedit\tupTIR\tdownTIR\tTIRedit\tnestedInside\tcontainNested" > elementFinder/penicillium.elements.ann.feat
join -t$'\t' -1 1 -2 1 <(sort -t$'\t' -k1,1 elementFinder/penicillium_elements_vs_YRsuperfams_besthits.txt | grep -P '_e|_s' | cut -f1,2) <(sort -t$'\t' -k1,1 elementFinder/penicillium.elements.temp.feat) | awk -F'\t' '{print}' >> elementFinder/penicillium.elements.ann.feat

mkdir pairViz
//start 13h20
starfish pair-viz -m all -t empty -T 2 -A nucmer -a ome2assembly.txt -b elementFinder/penicillium.elements.bed -S elementFinder/penicillium.elements.named.stats -o pairViz/

--------------FLANK--------------

mkdir flank_out

starfish flank   -a ome2assembly.txt   -b elementFinder/penicillium.insert.bed   -x macpha6   -o flank_out



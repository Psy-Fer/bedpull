Here are a few things I tried along the way to try and extract subsequences from the hg002 assembly

### commands to get regions from hg002

#### split the hg002 genome into maternal and paternal genomes
grep "_MATERNAL" hg002v1.1.fasta.fai | cut -f1 | xargs samtools faidx hg002v1.1.fasta > hg002.maternal_genome.fasta

grep "_PATERNAL" hg002v1.1.fasta.fai | cut -f1 | xargs samtools faidx hg002v1.1.fasta > hg002.paternal_genome.fasta

#### get genome sizes from hs1

cut -f1,2 hs1.fa.fai > hs1.genome.sizes

#### extract 50k left flanks

bedtools flank -i STRchive-disease-loci.T2T-chm13.TRGT.bed -g hs1.genome.sizes -l 50000 -r 0 > hs1.left_flanks.bed

#### extract 50k right flanks
bedtools flank -i STRchive-disease-loci.T2T-chm13.TRGT.bed -g hs1.genome.sizes -l 0 -r 50000 > hs1.right_flanks.bed

#### merge flanks together
cat hs1.left_flanks.bed hs1.right_flanks.bed > hs1.all_flanks.bed

#### get flank sequences
bedtools getfasta -fi hs1.fa -bed hs1.all_flanks.bed -fo hs1.flanks.fasta

#### map the sequences to the hg002 paternal/maternal genome
minimap2 -t 20 -ax asm5 hg002.paternal_genome.fasta hs1.flanks.fasta > flank_alignment_paternal.sam

#### Get the coords of the alignments

samtools view -bS -F 2304m flank_alignment_paternal.sam | bedtools bamtobed -i stdin | awk '{print $4, $1, $2, $3}' > flank_coords_paternal.txt

#### after extracting the coords (by hand, TODO: write a script) extract the sequence from the hg002 genome for the regions

bedtools getfasta -fi hg002.paternal_genome.fasta -bed paternal_cuts.bed -fo paternal_sequences.fasta -name

#### analyse the sequences for repeats to get baseline



# Whole Genome alignment method of getting co-ordinates

## HG002 aligned to hs1

### minimap2 (works waaaay better)

#### Split diploid genome 
grep "_MATERNAL" hg002v1.1.fasta.fai | cut -f1 | xargs samtools faidx hg002v1.1.fasta > hg002.maternal_genome.fasta
grep "_PATERNAL" hg002v1.1.fasta.fai | cut -f1 | xargs samtools faidx hg002v1.1.fasta > hg002.paternal_genome.fasta

#### Align with minimap2 (paf output)
minimap2 -cx asm5 --cs=long -t 16 hs1.fa hg002.maternal_genome.fasta > hg002mat_to_hs1.paf
minimap2 -cx asm5 --cs=long -t 16 hs1.fa hg002.paternal_genome.fasta > hg002pat_to_hs1.paf

#### Convert PAF to chain format (maf-convert is from LAST)
paftools.js view -f maf hg002mat_to_hs1.paf | maf-convert chain - > hg002mat_to_hs1.chain
paftools.js view -f maf hg002pat_to_hs1.paf | maf-convert chain - > hg002pat_to_hs1.chain

#### Lift coordinates
liftOver hs1_regions.bed hg002mat_to_hs1.chain hg002_maternal_regions.bed unmapped_mat.bed
liftOver hs1_regions.bed hg002pat_to_hs1.chain hg002_paternal_regions.bed unmapped_pat.bed

#### Results
  55 hg002_maternal_regions.bed
  52 hg002_paternal_regions.bed

### LAST

#### Split diploid genome
grep "_MATERNAL" hg002v1.1.fasta.fai | cut -f1 | xargs samtools faidx hg002v1.1.fasta > hg002.maternal_genome.fasta
grep "_PATERNAL" hg002v1.1.fasta.fai | cut -f1 | xargs samtools faidx hg002v1.1.fasta > hg002.paternal_genome.fasta

#### Build hs1 database
lastdb -P 16 -uNEAR hs1_db hs1.fa

#### Align and create chain files
lastal -P 16 -C2 hs1_db hg002.maternal_genome.fasta | last-split | maf-convert chain - > hg002mat_to_hs1.chain
lastal -P 16 -C2 hs1_db hg002.paternal_genome.fasta | last-split | maf-convert chain - > hg002pat_to_hs1.chain

#### Lift coordinates
liftOver hs1_regions.bed hg002mat_to_hs1.chain hg002_maternal_regions.bed unmapped_mat.bed
liftOver hs1_regions.bed hg002pat_to_hs1.chain hg002_paternal_regions.bed unmapped_pat.bed

#### Results
  15 hg002_maternal_regions.bed
  16 hg002_paternal_regions.bed




## got mad

wrote bedpull
got correct extractions
#!/bin/sh

#This a simplified script that shows the essential steps of the analysis


all_sets="set_A set_B set_C"
samples_A="ERR2224046 ERR2224047 ERR2224048 ERR2224049 ERR2224050 ERR2224051"
samples_B="SRR6960207 SRR6960209 SRR6960211 SRR6960213 SRR6960215 SRR6960217 SRR6960219 SRR6960221 SRR6960223 SRR6960225 SRR6960208 SRR6960210 SRR6960212 SRR6960214 SRR6960216 SRR6960218 SRR6960220 SRR6960222 SRR6960224"
samples_C="SRR3197701 SRR3197704 SRR3197707 SRR3197710 SRR3197700 SRR3197703 SRR3197706 SRR3197709 SRR3197702 SRR3197705 SRR3197708 SRR3197711"
all_samples="ERR2224046 ERR2224047 ERR2224048 ERR2224049 ERR2224050 ERR2224051 SRR6960207 SRR6960209 SRR6960211 SRR6960213 SRR6960215 SRR6960217 SRR6960219 SRR6960221 SRR6960223 SRR6960225 SRR6960208 SRR6960210 SRR6960212 SRR6960214 SRR6960216 SRR6960218 SRR6960220 SRR6960222 SRR6960224 SRR3197701 SRR3197704 SRR3197707 SRR3197710 SRR3197700 SRR3197703 SRR3197706 SRR3197709 SRR3197702 SRR3197705 SRR3197708 SRR3197711"




#for dataset C, the paired were not called R1 and R2 but just 1 and 2
for i in ${samples_C}
do
mv $work_dir/01-Fastq/set_C/${i}_1.fastq.gz $work_dir/01-Fastq/set_C/${i}_R1.fastq.gz
mv $work_dir/01-Fastq/set_C/${i}_2.fastq.gz $work_dir/01-Fastq/set_C/${i}_R2.fastq.gz
done


#======== Alignment

mkdir -p $work_dir/02-Alignment
mkdir -p $work_dir/02-Alignment/sam
mkdir -p $work_dir/02-Alignment/sam/set_A
mkdir -p $work_dir/02-Alignment/sam/set_B
mkdir -p $work_dir/02-Alignment/sam/set_C


for i in ${samples_A}
do
srun hisat2 -q --summary-file $work_dir/02-Alignment/sam/set_A/summary_${i}.txt -x /shared/projects/minomics/pgrognet/GenomePodo/index_hisat2/podo -S $work_dir/02-Alignment/sam/set_A/${i}.sam -1 #$work_dir/01-Fastq/set_A/${i}_R1.fastq.gz -2 $work_dir/01-Fastq/set_A/${i}_R2.fastq.gz
done


for i in ${samples_C}
do
srun hisat2 -q --summary-file $work_dir/02-Alignment/sam/set_C/summary_${i}.txt -x /shared/projects/minomics/pgrognet/GenomePodo/index_hisat2/podo -S $work_dir/02-Alignment/sam/set_C/${i}.sam -1 #$work_dir/01-Fastq/set_C/${i}_R1.fastq.gz -2 $work_dir/01-Fastq/set_C/${i}_R2.fastq.gz
done

for i in ${samples_B}
do
srun hisat2 -q --summary-file $work_dir/02-Alignment/sam/set_B/summary_${i}.txt -x /shared/projects/minomics/pgrognet/GenomePodo/index_hisat2/podo -S $work_dir/02-Alignment/sam/set_B/${i}.sam -U #$work_dir/01-Fastq/set_B/${i}.fastq.gz
done

echo "========= ALL ALIGNMENTS DONE ==========="


#make bam files

mkdir -p $work_dir/02-Alignment/bam/

for i in ${samples_A}
do
srun samtools sort $work_dir/02-Alignment/sam/set_A/${i}.sam | samtools view -b -o $work_dir/02-Alignment/bam/${i}.bam
srun samtools index $work_dir/02-Alignment/bam/${i}.bam
done

for i in ${samples_B}
do
srun samtools sort $work_dir/02-Alignment/sam/set_B/${i}.sam | samtools view -b -o $work_dir/02-Alignment/bam/${i}.bam
srun samtools index $work_dir/02-Alignment/bam/${i}.bam
done

for i in ${samples_C}
do
srun samtools sort $work_dir/02-Alignment/sam/set_C/${i}.sam | samtools view -b -o $work_dir/02-Alignment/bam/${i}.bam
srun samtools index $work_dir/02-Alignment/bam/${i}.bam
done


echo "===================== BAM files DONE ============="

#counting reads


for i in ${samples_A}
do
echo ${i} > $work_dir/01-Fastq/count_fastq_${i}.txt
expr $(cat $work_dir/01-Fastq/set_A/${i}_R1.fastq | wc -l) / 4 >> $work_dir/01-Fastq/count_fastq_${i}.txt
done

for i in ${samples_B}
do
echo ${i} > $work_dir/01-Fastq/count_fastq_${i}.txt
expr $(cat $work_dir/01-Fastq/set_B/${i}.fastq | wc -l) / 4 >> $work_dir/01-Fastq/count_fastq_${i}.txt
done

for i in ${samples_C}
do
echo ${i} > $work_dir/01-Fastq/count_fastq_${i}.txt
expr $(cat $work_dir/01-Fastq/set_C/${i}_R1.fastq | wc -l) / 4 >> $work_dir/01-Fastq/count_fastq_${i}.txt
done

echo "===================== count reads DONE ============="

#Conversion to bedgraph/bigwig for data visualization

mkdir -p $work_dir/03-visualization

for i in $all_samples
do
srun bamCoverage --bam $work_dir/02-Alignment/bam/${i}.bam -o $work_dir/02-Alignment/bam/${i}.bigwig -of bigwig --binSize 10 --normalizeUsing BPM
srun bamCoverage --bam $work_dir/02-Alignment/bam/${i}.bam -o $work_dir/02-Alignment/bam/${i}.bedgraph -of bedgraph --binSize 10 --normalizeUsing BPM
done
mv $work_dir/02-Alignment/bam/*.bigwig $work_dir/03-visualization
mv $work_dir/02-Alignment/bam/*.bedgraph $work_dir/03-visualization

echo "===================== BEDgraph/bigwig files DONE ============="




#StringTie

mkdir -p $work_dir/04-StringTie

for i in ${samples_A}
do
srun stringtie $work_dir/02-Alignment/bam/${i}.bam -o $work_dir/04-StringTie/${i}.gtf -c 10 -g 5
done

for i in ${samples_C}
do
srun stringtie $work_dir/02-Alignment/bam/${i}.bam -o $work_dir/04-StringTie/${i}.gtf -c 10 -g 5
done

for i in ${samples_B}
do
srun stringtie $work_dir/02-Alignment/bam/${i}.bam -o $work_dir/04-StringTie/${i}.gtf --rf -c 10 -g 5
done

echo "===================== stringtie GTF files DONE ============="

#StringTie merge

mkdir -p $work_dir/04-StringTie/merge

srun stringtie --merge $work_dir/04-StringTie/*.gtf -o $work_dir/04-StringTie/merge/merge_transcripts.gtf

echo "===================== stringtie merge DONE ============="





# HTseq count

mkdir -p $work_dir/05-htseqcount/

for i in ${samples_A}
do
srun htseq-count --stranded=no --type='gene' --idattr='gene_id' --order=name --format=bam -m intersection-nonempty $work_dir/02-Alignment/bam/${i}.bam $work_dir/Annotation_Podospora_Mat+_complet.gtf > $work_dir/05-htseqcount/count-${i}.txt
done

for i in ${samples_C}
do
srun htseq-count --stranded=no --type='gene' --idattr='gene_id' --order=name --format=bam -m intersection-nonempty $work_dir/02-Alignment/bam/${i}.bam $work_dir/Annotation_Podospora_Mat+_complet.gtf > $work_dir/05-htseqcount/count-${i}.txt
done

for i in ${samples_B}
do
srun htseq-count --stranded=reverse --type='gene' --idattr='gene_id' --order=name --format=bam -m intersection-nonempty $work_dir/02-Alignment/bam/${i}.bam $work_dir/Annotation_Podospora_Mat+_complet.gtf > $work_dir/05-htseqcount/count-${i}.txt
done

#echo "===================== HTseq count files DONE =================="





#TopHat2

mkdir -p $work_dir/07-TopHat2/

tophat --min-intron-length 30 --max-multihits 5 /shared/projects/minomics/pgrognet/GenomePodo/index $work_dir/01-Fastq/set_C/
tophat --min-intron-length 30 --max-multihits 5 --segment-length 21 /shared/projects/minomics/pgrognet/GenomePodo/index $work_dir/01-Fastq/set_A
tophat --min-intron-length 30 --max-multihits 5 --library-type fr-firstrand /shared/projects/minomics/pgrognet/GenomePodo/index $work_dir/01-Fastq/set_B
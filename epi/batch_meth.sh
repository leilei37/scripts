prefetch=/public/home/leimengyu/atac/sratoolkit.3.0.0-ubuntu64/bin/prefetch
fastq_dump=/public/home/leimengyu/epi/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump
fastqc=/public/home/muvunyi/miniconda3/bin/fastqc
multiqc=/public/home/muvunyi/miniconda3/bin/multiqc
bowtie2=/public/software/bin/bowtie2
samtools=/public/software/bin/samtools
trimmomatic=~/anaconda3/bin/trimmomatic
wkdir='/public/home/leimengyu/epi/meth'
picard=~/anaconda3/bin/picard

for line in `cat /public/home/leimengyu/epi/meth/sample.txt`;do
#data list  collection: ncbi search ‘sra B73 epigenome’
nohup $prefetch SRR &
${fastq_dump} --split-3 --gzip ${wkdir}/${line}/${line}.sra
#trim
$trimmomatic PE     -threads 16     -phred33     -trimlog trim.log     /public/home/leimengyu/${line}_1.fastq.gz /public/home/leimengyu/${line}_2.fastq.gz     ${wkdir}/trim/${line}_R1.fq.gz ${wkdir}/trim/${line}_unpaired_R1.fq.gz ${wkdir}/trim/${line}_R2.fq.gz ${wkdir}/trim/${line}_unpaired_R2.fq.gz     ILLUMINACLIP:/public/home/leimengyu/epi/atac/NexteraPE-PE.fa:2:30:10:8:true     SLIDINGWINDOW:4:15     LEADING:3     TRAILING:3 
#fastqc
$fastqc --outdir ${wkdir}/fastqc/ --threads 16 ${wkdir}/trim/${line}_R1.fq.gz
$fastqc --outdir ${wkdir}/fastqc/ --threads 16 ${wkdir}/trim/${line}_R2.fq.gz
#align
bismark --parallel 16 --non_directional -o $wkdir/trim/${line}bismark_file_name_trimmed.n5 -N 1 --genome /public/home/leimengyu/data/methgeno/ -1 $wkdir/trim/${line}_R1.fq.gz -2 $wkdir/trim/${line}_R2.fq.gz
#remove deduplicate，-p for pair-end data
deduplicate_bismark -p $wkdir/trim/${line}bismark_file_name_trimmed.n5/${line}_R1_bismark_bt2_pe.bam --output_dir $wkdir/trim/${line}bismark_file_name_trimmed.n5/ --samtools_path /public/home/leimengyu/anaconda3/envs/epi/bin
#extract
bismark_methylation_extractor --gzip --bedGraph $wkdir/trim/${line}bismark_file_name_trimmed.n5/${line}_R1_bismark_bt2_pe.deduplicated.bam -o $wkdir/trim/${line}bismark_file_name_trimmed.n5/
#bismark_methylation_extractor -p --comprehensive --no_overlap --bedGraph --counts --buffer_size 20G --report --cytosine_report --samtools_path /cm/shared/opt/anaconda2/bin/ --genome_folder /data/ref_fastq_2/ /result/S.deduplicated.sam.gz -o /result
done

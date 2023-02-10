prefetch=/public/home/leimengyu/atac/sratoolkit.3.0.0-ubuntu64/bin/prefetch
${fastq-dump}=/public/home/leimengyu/epi/sratoolkit.3.0.0-ubuntu64/bin/fastq-dump
fastqc=/public/home/muvunyi/miniconda3/bin/fastqc
multiqc=/public/home/muvunyi/miniconda3/bin/multiqc
bowtie2=/public/software/bin/bowtie2
samtools=/public/software/bin/samtools
trimmomatic=~/anaconda3/bin/trimmomatic
wkdir='/public/home/leimengyu/epi/atac'
picard=~/anaconda3/bin/picard
gatk=~/anaconda3/bin/gatk
macs2=~/anaconda3/envs/python3.6/bin/macs2

for line in `cat /public/home/leimengyu/epi/atac/sample.txt`;do
#data list  collection: ncbi search ‘sra B73 epigenome’
#nohup $prefetch SRR &
${fastq-dump} --split-3 --gzip ${wkdir}/${line}/${line}.sra
#trim
$trimmomatic PE     -threads 16     -phred33     -trimlog trim.log     ${wkdir}/${line}/${line}_1.fastq.gz ${wkdir}/${line}/${line}_2.fastq.gz     ${wkdir}/trim/${line}_R1.fq.gz ${wkdir}/trim/${line}_unpaired_R1.fq.gz ${wkdir}/trim/${line}_R2.fq.gz ${wkdir}/trim/${line}_unpaired_R2.fq.gz     ILLUMINACLIP:${wkdir}/NexteraPE-PE.fa:2:30:10:8:true     SLIDINGWINDOW:5:20     LEADING:3     TRAILING:3     MINLEN:30
#fastqc
$fastqc --outdir ${wkdir}/fastqc/ --threads 16 ${wkdir}/trim/${line}_R1.fq.gz
$fastqc --outdir ${wkdir}/fastqc/ --threads 16 ${wkdir}/trim/${line}_R2.fq.gz
#align
$bowtie2  -p 16 --very-sensitive  -X 2000  -x /public/home/leimengyu/data/maizeindex  -1 ${wkdir}/trim/${line}_R1.fq.gz  -2 ${wkdir}/trim/${line}_R2.fq.gz -S ${wkdir}/bowtie2/${line}.sam
#summary for align
$samtools view  -bS  ${wkdir}/bowtie2/${line}.sam > ${wkdir}/bowtie2/${line}.1bam
$samtools sort ${wkdir}/bowtie2/${line}.1bam ${wkdir}/bowtie2/${line}.bam
$samtools index ${wkdir}/bowtie2/${line}.bam
$samtools flagstat ${wkdir}/bowtie2/${line}.bam > ${wkdir}/bowtie2/${line}.stat
$samtools idxstats ${wkdir}/bowtie2/${line}.bam > ${wkdir}/bowtie2/${line}_idxstats.txt
$gatk CollectInsertSizeMetrics -H ${wkdir}/bowtie2/${line}_InsertSize.pdf -I ${wkdir}/bowtie2/${line}.bam -O ${wkdir}/bowtie2/${line}_InsertSize.txt
#remove duplicates
$picard MarkDuplicates REMOVE_DUPLICATES=true I=${wkdir}/bowtie2/${line}.bam O=${wkdir}/bowtie2/${line}markdup.bam M=${wkdir}/bowtie2/${line}markdup.metrc.csvp
#filter low quality aligned reads and MT
$samtools view -h ${wkdir}/bowtie2/${line}markdup.bam|grep -v Mt|$samtools view -bS -q 30 -F 1804 -f 0x2 - >${wkdir}/bowtie2/${line}filt.bam
#call peak
$macs2 callpeak -f BAMPE -t ${wkdir}/bowtie2/${line}filt.bam  -g 2.1e+9 --nomodel --shift -50  --extsize 100  -n ${line} -B -q 0.01 --keep-dup all --outdir ${wkdir}/macs3/
#Fraction of reads in peaks (FRiP)
$samtools view -c ${wkdir}/bowtie2/${line}filt.bam > ${wkdir}/bowtie2/${line}.count
$samtools view -c -L ${wkdir}/macs3/${line}.narrowPeak {wkdir}/bowtie2/${line}filt.bam > ${wkdir}/bowtie2/${line}.peakcount
done
idr --samples peak1 peak2 --peak-list merge.peak --plot

#make sample.txt in ${wkdir}
fastqc='/public/home/leimengyu/anaconda3/bin/fastqc'
multiqc='/public/home/muvunyi/miniconda3/bin/multiqc'
STAR='/public/home/leimengyu/anaconda3/bin/STAR'
qualimap='/public/home/leimengyu/anaconda3/bin/qualimap'
gffread='/public/home/leimengyu/anaconda3/bin/gffread'
featureCounts='/public/home/leimengyu/anaconda3/bin/featureCounts'
wkdir='/public/home/leimengyu/mengyun/'
sample='/public/home/leimengyu/mengyun/sample.txt'
sampledir='/public/home/leimengyu/ANNO_ANCBJ170219_PM-ANCBJ170219-13_2022-06-06_14-46-37_BH33FCDSX3/Cleandata'
#####sequence qc 
for line in `cat $sample`;do
$fastqc -t 16 -o $wkdir ${sampledir}/${line}/*.fq.gz
done
$multiqc ${wkdir}/*zip -o $wkdir
###make index(for once)
#$STAR --runMode genomeGenerate --genomeDir /public/home/leimengyu/riceg/osindex/ --genomeFastaFiles /public/home/leimengyu/riceg/all.con --sjdbGTFfile /public/home/leimengyu/riceg/all.gff3 --sjdbOverhang 150 --genomeSAindexNbases 10
#####align
for line in `cat $sample`;do
$STAR --runThreadN 40 --genomeDir /public/home/leimengyu/riceg/osindex --readFilesCommand zcat --readFilesIn ${sampledir}/${line}/${line}_R1.fq.gz ${sampledir}/${line}/${line}_R2.fq.gz --outFileNamePrefix ${wkdir}bam/$line --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 40 --quantMode TranscriptomeSAM GeneCounts 
done
###gff convert to gtf(for once)
#$gffread /public/home/leimengyu/riceg/all.gff3 -T -o /public/home/leimengyu/riceg/all.gtf
#####make gene expression matrix
$featureCounts -a /public/home/leimengyu/riceg/all.gtf -p -T 8 -o ${wkdir}matrix.txt ${wkdir}bam/*.sortedByCoord.out.bam
#####make mapfile for bamqc
ls ${wkdir}bam/*Aligned.sort* > ${wkdir}bamfile.txt
paste -d "\t" ${wkdir}sample.txt ${wkdir}bamfile.txt > ${wkdir}mapfile.txt
#####bamqc
$qualimap multi-bamqc -r -d ${wkdir}mapfile.txt -gff /public/home/leimengyu/riceg/all.gtf -outdir ${wkdir}multi_bamqc -outformat PDF:HTML
#####summary uniquemapgene multimapgene
#in bam file
#python3 finduniquemapmultimap.py

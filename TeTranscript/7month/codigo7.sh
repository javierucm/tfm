
#Primer paso: Descarga dropbox GTF TE
#wget https://ftp.ensembl.org/pub/release-112/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
#wget https://ftp.ensembl.org/pub/release-112/gtf/mus_musculus/Mus_musculus.GRCm39.112.gtf.gz
#gunzip Mus_musculus.GRCm39.dna.primary_assembly.fa.gz 

#STAR --runMode genomeGenerate --genomeDir genomeTE --genomeFastaFiles Mus_musculus.GRCm39.dna.primary_assembly.fa --runThreadN 15




######################## Samples information and indentity ########################

#SRR17261749	GSM5737122: Young_Skin_non-dox_5287; Mus musculus; RNA-Seq
#SRR17261748    GSM5737121: Young_Skin_non-dox_5281; Mus musculus; RNA-Seq
#SRR17261747	GSM5737120: Young_Skin_non-dox_5275; Mus musculus; RNA-Seq
#SRR17261746	GSM5737119: Young_Skin_non-dox_5269; Mus musculus; RNA-Seq

#SRR17261674 	GSM5737047	Skin_4F_6545
#SRR17261675 	GSM5737048	Skin_4F_6546
#SRR17261676 	GSM5737049	Skin_4F_6547
#SRR17261677 	GSM5737050	Skin_4F_6548
#SRR17261678	GSM5737051	Skin_4F_6549

######################## Download samples ########################

#fastq-dump --gzip --split-3 -v SRR17261749
#fastq-dump --gzip --split-3 -v SRR17261748
#fastq-dump --gzip --split-3 -v SRR17261747
#fastq-dump --gzip --split-3 -v SRR17261746
#fastq-dump --gzip --split-3 -v SRR17261674
#fastq-dump --gzip --split-3 -v SRR17261675
#fastq-dump --gzip --split-3 -v SRR17261676
#fastq-dump --gzip --split-3 -v SRR17261677
#fastq-dump --gzip --split-3 -v SRR17261678


for sample in SRR17261749 SRR17261748 SRR17261747 SRR17261746 SRR17261674 SRR17261675 SRR17261676 SRR17261677 SRR17261678
do
STAR --runMode alignReads --runThreadN 7 --genomeDir genomeTE/ --sjdbGTFfile Mus_musculus.GRCm39.112.gtf --outSAMtype BAM SortedByCoordinate --winAnchorMultimapNmax 100 --outFilterMultimapNmax 100 --readFilesCommand zcat --outFileNamePrefix ${sample}/${sample}_ --readFilesIn ${sample}.fastq.gz
done

#OJO STRANDEDDDDDDDDDDD

# Extraer los exones del archivo GTF y convertirlos a BED
#awk '$3 == "exon"' Mus_musculus.GRCm39.112.gtf | \
#awk 'BEGIN{OFS="\t"} {print $1, $4-1, $5, ".", ".", $7}' > genes.bed



#infer_experiment.py -r genes.bed  -i SRR17261749_Aligned.sortedByCoord.out.bam

#--reverse stranded

# nohup TEtranscripts --format BAM --stranded reverse -t SRR17261674/SRR17261674_Aligned.sortedByCoord.out.bam SRR17261675/SRR17261675__Aligned.sortedByCo
ord.out.bam SRR17261676/SRR17261676_Aligned.sortedByCoord.out.bam SRR17261677/SRR17261677_Aligned.sortedByCoord.out.bam SRR17261678/SRR17261678_Aligned.sortedByCoord.out.bam -c SRR17261746/SRR17261746_Aligned.sortedByCoord.
out.bam  SRR17261747/SRR17261747_Aligned.sortedByCoord.out.bam SRR17261748/SRR17261748_Aligned.sortedByCoord.out
.bam SRR17261749/SRR17261749_Aligned.sortedByCoord.out.bam  --GTF Mus_musculus.GRCm39.112.gtf --TE GRCm39_Ensemb
l_rmsk_TE.gtf --mode multi --sortByPos --project TFM_JLR_Long7MonthYoungvsOldTreated



#Los genes significativamente diferenciados me salen distintos, en número, con una gran diferencia, cuando hago el Deseq sobre los 10, que cuando lo hago sobre las 7 replicas. 
# hay un número grande.

#hago el intersect o la union?




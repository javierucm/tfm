######################## Samples information and indentity ########################

#SRR17261749	GSM5737122: Young_Skin_non-dox_5287; Mus musculus; RNA-Seq
#SRR17261748    GSM5737121: Young_Skin_non-dox_5281; Mus musculus; RNA-Seq
#SRR17261747	GSM5737120: Young_Skin_non-dox_5275; Mus musculus; RNA-Seq
#SRR17261746	GSM5737119: Young_Skin_non-dox_5269; Mus musculus; RNA-Seq

#SRR17261723	GSM5737096: Old_Skin_non-dox_5245; Mus musculus; RNA-Seq
#SRR17261722	GSM5737095: Old_Skin_non-dox_5239; Mus musculus; RNA-Seq
#SRR17261721	GSM5737094: Old_Skin_non-dox_5233; Mus musculus; RNA-Seq

#SRR17261720	GSM5737093: Old_Skin_dox-treated_5263; Mus musculus; RNA-Seq
#SRR17261719	GSM5737092: Old_Skin_dox-treated_5257; Mus musculus; RNA-Seq
#SRR17261718	GSM5737091: Old_Skin_dox-treated_5251; Mus musculus; RNA-Seq
######################## Download samples ########################

#fastq-dump --gzip --split-3 -v SRR17261749
#fastq-dump --gzip --split-3 -v SRR17261748
#fastq-dump --gzip --split-3 -v SRR17261747
#fastq-dump --gzip --split-3 -v SRR17261746



#wget https://ftp.ensembl.org/pub/release-112/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz

#STAR --runMode genomeGenerate --genomeDir genomeTE --genomeFastaFiles Mus_musculus.GRCm39.dna.primary_assembly.fa --runThreadN 15

#Descarga dropbox GTF


#wget https://ftp.ensembl.org/pub/release-112/gtf/mus_musculus/Mus_musculus.GRCm39.112.gtf.gz




#for sample in SRR17261723 SRR17261722 SRR17261721
for sample in  SRR17261719
do
STAR --runMode alignReads --runThreadN 7 --genomeDir genomeTE/ --sjdbGTFfile Mus_musculus.GRCm39.112.gtf --outSAMtype BAM SortedByCoordinate --winAnchorMultimapNmax 100 --outFilterMultimapNmax 100 --readFilesCommand zcat --outFileNamePrefix ${sample}/${sample}_ --readFilesIn ${sample}.fastq.gz
done

#OJO STRANDEDDDDDDDDDDD

# Extraer los exones del archivo GTF y convertirlos a BED
awk '$3 == "exon"' Mus_musculus.GRCm39.112.gtf | \
awk 'BEGIN{OFS="\t"} {print $1, $4-1, $5, ".", ".", $7}' > genes.bed



infer_experiment.py -r genes.bed  -i SRR17261749_Aligned.sortedByCoord.out.bam

--reverse stranded

# TEtranscripts --format BAM --stranded reverse -t SRR17261723/SRR17261723_Aligned.sortedByCoord.out.bam SRR17261722/SRR17261722_Aligned.sortedByCoord.out.bam SRR17261721/SRR17261721_Aligned.sortedByCoord.out.bam -c SRR17261749_Aligned.sortedByCoord.out.bam  SRR17261748/SRR17261748_Aligned.sortedByCoord.out.bam SRR17261747/SRR17261747_Aligned.sortedByCoord.out.bam SRR17261746/SRR17261746_Aligned.sortedByCoord.out.bam  --GTF Mus_musculus.GRCm39.112.gtf --TE GRCm39_Ensembl_rmsk_TE.gtf --mode multi --sortByPos --project TFM_Javier_Lazaro



#Los genes significativamente diferenciados me salen distintos, en número, con una gran diferencia, cuando hago el Deseq sobre los 10, que cuando lo hago sobre las 7 replicas. 
# hay un número grande.

#hago el intersect o la union?




# **Genomic analysis of 2023 Summer AZ field data**

**Author: Li'ang Yu, Andrew Nelson** \
**Date: Dec 13th, 2023**

- [**Genomic analysis of 2023 Summer AZ field cotton transcriptome data**](#genomic-analysis-of-2023-summer-az-field-cotton-transcriptome-data)
  - [ATAC-seq data trail experiemnts](#atac-seq-data-trail-experiemnts)
    - [Trim reads using trimGalore](#trim-reads-using-trimgalore)
    - [Check trimmed fastq reads numbers](#check-trimmed-fastq-reads-numbers)
    - [Mapping reads to the genome (With two different mapping strategy)](#mapping-reads-to-the-genome-with-two-different-mapping-strategy)
      - [A and D-subgenome](#a-and-d-subgenome)
      - [A and D (whole genome combined)](#a-and-d-whole-genome-combined)
      - [Entire cotton genome](#entire-cotton-genome)
  - [Transcriptome data processing](#transcriptome-data-processing)
    - [Trim raw reads](#trim-raw-reads)


## ATAC-seq data trail experiemnts
### Trim reads using trimGalore
```bash
#! /bin/bash
raw_reads_dir="/mnt/Dominique/01_Cotton_2023/01_ATAC-trail"
trim_reads_dir="/mnt/Dominique/01_Cotton_2023/01_ATAC-trail/02_ATAC-trim"

/home/xzhang/Software/TrimGalore/trim_galore \
  --phred33 \
  --paired \
  --fastqc \
  --cores 16 \
  --output_dir /mnt/Dominique/01_Cotton_2023/01_ATAC-trail/02_ATAC-trim \
  --path_to_cutadapt /home/liangyu/anaconda3/envs/cutadaptenv/bin/cutadapt \
  ${raw_reads_dir}/Yinxin_cotton2-1_S9_L003_R1_001.fastq.gz ${raw_reads_dir}/Yinxin_cotton2-1_S9_L003_R2_001.fastq.gz

/home/xzhang/Software/TrimGalore/trim_galore \
  --phred33 \
  --paired \
  --fastqc \
  --cores 16 \
  --output_dir /mnt/Dominique/01_Cotton_2023/01_ATAC-trail/02_ATAC-trim \
  --path_to_cutadapt /home/liangyu/anaconda3/envs/cutadaptenv/bin/cutadapt \
  ${raw_reads_dir}/Yinxin_cotton2-2_S10_L003_R1_001.fastq.gz ${raw_reads_dir}/Yinxin_cotton2-2_S10_L003_R2_001.fastq.gz
```

### Check trimmed fastq reads numbers
```bash
echo $(zcat Coker_310_WL_A_1_val_1.fq.gz | wc -l)/4|bc 
echo $(zcat Coker_310_WL_A_1_val_1.fq.gz | wc -l)/4|bc 
echo $(zcat Coker_310_WL_A_1_val_1.fq.gz | wc -l)/4|bc 
echo $(zcat Coker_310_WL_A_1_val_1.fq.gz | wc -l)/4|bc 
```

Output and sample name
```
echo $(zcat Coker_310_WL_A_1_val_1.fq.gz | wc -l)/4|bc    34,440,608
echo $(zcat oker_310_WL_A_forward_paired.fq.gz | wc -l)/4|bc     7,473,506
```

### Mapping reads to the genome (With two different mapping strategy)
#### A and D-subgenome 
```bash
mapped_reads_dir="/mnt/Dominique/01_Cotton_2023/02_ATAC-bam"

#!/bin/bash
Ref_dir="/mnt/Knives/1_data/9_Cotton_omics/0_data"
reads_dir="/mnt/Dominique/01_Cotton_2023/01_ATAC-trail/02_ATAC-trim" 
workdir="/mnt/Dominique/01_Cotton_2023/02_ATAC-bam/01_Subgenome" 

#Build index and dict files for genome
samtools faidx $Ref_dir/Ghirsutum_527_A_sub.fa
bwa index $Ref_dir/Ghirsutum_527_A_sub.fa 
picard CreateSequenceDictionary -R $Ref_dir/Ghirsutum_527_A_sub.fa -O $Ref_dir/Ghirsutum_527_A_sub.dict

#for the A-subgenome
for file in $(cat ${reads_dir}/sample.list);
do

  ## Map reads by bwa piepline
  bwa mem -M -R "@RG\tID:$file\tSM:$file\tPL:ILLUMINA" \
    -t 14 -B 40 \
    $Ref_dir/Ghirsutum_527_A_sub.fa \
    ${reads_dir}/${file}_1_val_1.fq.gz \
    ${reads_dir}/${file}_2_val_2.fq.gz | samtools view -@10 -bS - -o ${workdir}/${file}.subA.reorder.bam 

  ## uniquely mapped reads 
  sambamba view -t 15 -h -f bam \
    -F "mapping_quality >= 1 and not (unmapped or secondary_alignment) and not ([XA] != null or [SA] != null)" \
    ${workdir}/${file}.subA.reorder.bam \
    -o ${workdir}/${file}.subA.uniqmapped.bam

 ### Generate the stats for two bam files
 samtools flagstat ${workdir}/${file}.subA.reorder.bam > ${workdir}/${i}subA.reorder.stat.txt
 samtools flagstat ${workdir}/${file}.subA.uniqmapped.bam > ${workdir}/${i}.subA.uniqmapped.stat.txt
done

#for the D-subgenome
for file in $(cat ${reads_dir}/sample.list);
do

  ## Map reads by bwa piepline
  bwa mem -M -R "@RG\tID:$file\tSM:$file\tPL:ILLUMINA" \
    -t 14 -B 40 \
    $Ref_dir/Ghirsutum_527_D_sub.fa \
    ${reads_dir}/${file}_1_val_1.fq.gz \
    ${reads_dir}/${file}_2_val_2.fq.gz | samtools view -@10 -bS - -o ${workdir}/${file}.subD.reorder.bam 

  ## uniquely mapped reads 
  sambamba view -t 15 -h -f bam \
    -F "mapping_quality >= 1 and not (unmapped or secondary_alignment) and not ([XA] != null or [SA] != null)" \
    ${workdir}/${file}.subD.reorder.bam \
    -o ${workdir}/${file}.subD.uniqmapped.bam

 ### Generate the stats for two bam files
 samtools flagstat ${workdir}/${file}.subD.reorder.bam > ${workdir}/${i}.subD.reorder.stat.txt
 samtools flagstat ${workdir}/${file}.subD.uniqmapped.bam > ${workdir}/${i}.subD.uniqmapped.stat.txt
done
```

#### A and D (whole genome combined)
```bash
for file in $(cat ${reads_dir}/sample.list);
do

  ## Map reads by bwa piepline
  bwa mem -M -R "@RG\tID:$file\tSM:$file\tPL:ILLUMINA" \
    -t 14 -B 40 \
    $Ref_dir/Ghirsutum_527_v2.0.fa \
    ${reads_dir}/${file}_1_val_1.fq.gz \
    ${reads_dir}/${file}_2_val_2.fq.gz | samtools view -@10 -bS - -o ${workdir}/${file}.AD.reorder.bam 

  ## uniquely mapped reads 
  sambamba view -t 15 -h -f bam \
    -F "mapping_quality >= 1 and not (unmapped or secondary_alignment) and not ([XA] != null or [SA] != null)" \
    ${workdir}/${file}.AD.reorder.bam \
    -o ${workdir}/${file}.AD.uniqmapped.bam

 ### Generate the stats for two bam files
 samtools flagstat ${workdir}/${file}.AD.reorder.bam > ${workdir}/${i}.AD.reorder.stat.txt
 samtools flagstat ${workdir}/${file}.AD.uniqmapped.bam > ${workdir}/${i}.AD.uniqmapped.stat.txt
done
```

#### Entire cotton genome

## Transcriptome data processing
### Trim raw reads
```bash
RNA_reads_dir="/mnt/Dominique/01_Cotton_2023"
```


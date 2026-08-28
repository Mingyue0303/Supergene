### 1. isoseq3处理原始数据

```bash
# 1. skera 和 lima分样
skera split m84247_260327_192516_s4.hifi_reads.bcM0001.bam Mas8.adapter.fasta m84247_260327_192516_s4.hifi_reads.bcM0001.skera.bam
#2. 
lima -j 9 m84247_260327_192516_s4.hifi_reads.bcM0001.skera.bam Isoseq.V2.adapter m84247_260327_192516_s4.fl.bam --isoseq --peek-guess --overwrite-biosample-names

# 3. 去除adapter
isoseq3 refine --require-polya m84247_260327_192516_s4.fl.IsoSeqX_bc01_5p--IsoSeqX_3p.bam Isoseq.V2.adapter m84247_260327_192516_s4.fl.IsoSeqX_bc01_5p--IsoSeqX_3p.flnc.bam
...
# 4. mapping
pbmm2 align genome.fasta.mmi m84247_260327_192516_s4.fl.IsoSeqX_bc01_5p--IsoSeqX_3p.bam.flnc.bam m84247_260327_192516_s4.fl.IsoSeqX_bc01_5p--IsoSeqX_3p.flnc.bam.pbmm2.bam --preset ISOSEQ --sort
...
```

### 2. 处理BAM文件，对高覆盖度区间进行处理
```bash
1. 统计覆盖深度
module load SAMtools/1.22.1-GCC-14.2.0
samtools depth -a -Q 0 m84247_260327_192516_s4.fl.IsoSeqX_bc01_5p--IsoSeqX_3p.FLNC.bam.pbmm2.bam > coverage.txt
2. 筛选超过区域
module load BEDTools/2.30.0-GCC-10.3.0
awk '$3 > 5000 {print $1"\t"($2-1)"\t"$2}' coverage.txt > high_cov.bed
bedtools merge -i high_cov.bed -d 100 > high_cov_merged.bed
3. 提取高覆盖区间的 reads，随机下采样到 10%
module load SAMtools/1.22.1-GCC-14.2.0
module load BEDTools/2.30.0-GCC-10.3.0
samtools view -b -L high_cov_merged.bed m84247_260327_192516_s4.fl.IsoSeqX_bc01_5p--IsoSeqX_3p.FLNC.bam.pbmm2.bam \
  | samtools view -b -s 0.1 > downsampled_hotspot.bam

# 提取正常区间的 reads
samtools view -L high_cov_merged.bed m84247_260327_192516_s4.fl.IsoSeqX_bc12_5p--IsoSeqX_3p.FLNC.bam.pbmm2.bam | cut -f1 | sort -u > high.reads
samtools view m84247_260327_192516_s4.fl.IsoSeqX_bc12_5p--IsoSeqX_3p.FLNC.bam.pbmm2.bam | awk 'NR==FNR{a[$1];next} !($1 in a)' high.reads - > normal.sam
samtools view -H m84247_260327_192516_s4.fl.IsoSeqX_bc12_5p--IsoSeqX_3p.FLNC.bam.pbmm2.bam > header.sam
cat header.sam normal.sam | samtools view -b > normal.bam

# 合并
samtools merge -f merged_clean.bam normal.bam downsampled_hotspot.bam
samtools sort -o merged_clean_sorted.bam merged_clean.bam
samtools index merged_clean_sorted.bam

```
### 3.TAMA collapse
```bash

# 1. tama collapse 可根据文献设置low 和 high模式
source activate tama #需要python2.7环境
# 分别对每个样品跑tama collapse，下面以bc01样品为例
python /scratch/leuven/382/vsc38210/03container/tama/tama_collapse.py \
        -d merge_dup \
        -x no_cap \
        -a 100 \
        -z 100 \
        -sj sj_priority \
        -lde 1 \
        -sjt 20 \
        -b BAM \
        -f genome.fasta \
        -s merged_clean_sorted.bam \
        -p bc01 1>tama.log 2>&1
...
# 3. tama merge 因为我的isoform数量太多所以我设置了更为宽松的合并流程
source activate tama
python /scratch/leuven/382/vsc38210/03container/tama/tama_merge.py -f file.list -p merge -d merge_dup -m 20 -a 100 -z 100

# cat file.list
bc01.bed        no_cap  2,1,1   bc01
bc02.bed        no_cap  2,1,1   bc02
bc08.bed        no_cap  2,1,1   bc08
bc09.bed        no_cap  2,1,1   bc09
bc10.bed        no_cap  2,1,1   bc10
bc11.bed        no_cap  2,1,1   bc11
bc12.bed        no_cap  2,1,1   bc12
...
# 结果文件
merge.bed, merge_gene_report.txt, merge_merge.txt, merge_trans_report.txt
# 4. 计算read support
awk '{print $1"\t"$1"_trans_read.bed\ttrans_read"}' samples.txt > trans_read.fofn

cat trans_read.fofn
bc01    bc01_trans_read.bed     trans_read
bc02    bc02_trans_read.bed     trans_read
bc08    bc08_trans_read.bed     trans_read
bc09    bc09_trans_read.bed     trans_read
bc10    bc10_trans_read.bed     trans_read
bc11    bc11_trans_read.bed     trans_read
bc12    bc12_trans_read.bed     trans_read

python /gpfs1/scratch/382/vsc38210/03container/tama/tama_go/read_support/tama_read_support_levels.py -f trans_read.fofn -m merge_merge.txt -o merge

# 5. tama filter(默认参数）
python /scratch/leuven/382/vsc38210/03container/tama/tama_go/filter_transcript_models/tama_remove_single_read_models_levels.py \
        -b merge.bed \
        -r merge_read_support.txt \
        -o tama.filter.2.trans \
        -l transcript \
        -k remove_multi \
        -s 1 \
        -n 2
```
### 4. sqanti3 过滤
```
conda activate sqanti3
export LD_LIBRARY_PATH=/data/leuven/382/vsc38210/miniconda3/envs/sqanti3/lib:$LD_LIBRARY_PATH
python sqanti3_qc.py --isoforms tama22.trans.gtf \
        --refGTF isoquant.gtf \
        --refFasta genome.fasta \
        --short_reads short_reads.fofn \
        --fl_count fl_counts.tsv \
        --saturation \
        --cpus 10

python sqanti3_filter.py rules \
        --sqanti_class isoforms_classification.txt\
        --filter_gtf isoforms_corrected.gtf \
        --filter_isoforms isoforms_corrected.fasta \
        --filter_faa isoforms_corrected.faa
        -d . \
        -o tama2.sq\
        -j filter_default.json
# result file: tama2.sq.filtered.gtf
# HSS result file: hunch.tama2.sq.gtf
```
### 5. pasa 预测完整ORF
```
# Use the SQANTI3 output files as the PASA input files.
singularity exec -e -B /lustre1/scratch/382/vsc38210:/lustre1/scratch/382/vsc38210 --pwd $BASE $IMG /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl \
     -c $BASE/alignAssembly.config -C -R -g $BASE/genome.fasta  \
     -t $BASE/transcript.fasta \
     --ALIGNERS minimap2 --CPU 10  1>$BASE/pasa.log 2>&1

# 使用transdecoder获得读码框
singularity exec -e -B /lustre1/scratch/382/vsc38210:/lustre1/scratch/382/vsc38210 --pwd $BASE $IMG /usr/local/src/PASApipeline/scripts/pasa_asmbls_to_training_set.dbi \
        --pasa_transcripts_fasta $BASE/Ogib-pasa.sqlite.assemblies.fasta \
        --pasa_transcripts_gff3 $BASE/Ogib-pasa.sqlite.pasa_assemblies.gff3 \
        >$BASE/pasa2train.log 2>&1

# 提取读码框完整的结果
grep 'complete' Ogib-pasa.sqlite.assemblies.fasta.transdecoder.cds|sed 's/>\(\S\+\).*/\1/' > pasa.cds.complete.lst
grep -f pasa.cds.complete.lst Ogib-pasa.sqlite.assemblies.fasta.transdecoder.genome.gff3 > pasa.cds.complete.gff3
# HSS result file: hunch.tama2.pasa.gff3
```

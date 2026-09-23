### 1.salmon
```bash
# 1.运行salmon
conda activate salmon
salmon index -t braker4.transcripts.fa -i Ogib -p 16
salmon quant -i Ogib -l A -1 aaa_1.fq.gz -2 aaa_2.fq.gz -p 6 -o aaa

# 2.整理结果

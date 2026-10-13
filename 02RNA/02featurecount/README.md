### featurec caculates RNA-seq expression reads
```bash
# 1. running
conda activate RNA
featureCounts -p -T 5 -t exon -g gene_id -M \
    -a isoquant.gtf \
    -o  result/aaa.txt \
    ./bam/aaa

# 2.manage results
cat ***.bam.txt |cut -f1,7 > ***.bam.txt.count
paste all.count.file > Ogib.ftc.iso.all.csv
```

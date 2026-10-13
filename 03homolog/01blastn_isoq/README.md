### func.gene __ isoquant.transcrit.fa
```bash
# 1.search
module load BLAST+/2.16.0-gompi-2024a
makeblastdb \
    -in isoquant.trans.fa \
    -dbtype nucl \
    -parse_seqids \
    -out isoquant_db
blastn \
    -task blastn \
    -query all.hunch.func.trans.fa \
    -db isoquant_db \
    -out blastn_isoquant.tsv \
    -outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen" \
    -evalue 1e-10 \
    -num_threads 16
# 2.filter
python final.filterlen.uniq.py candidate_transcripts.txt isoquant_tx2gene.tsv blastn_isoquant.txt --min-qcov 10
result: blast_best_per_target_gene.tsv blast_pass_qcov.tsv
```

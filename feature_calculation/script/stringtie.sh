STAR --outFilterType BySJout --runThreadN 8 --genomeDir hg38_mapIndex --readFilesIn SRR500121.fastq --outFileNamePrefix SRR500121 --outSAMtype BAM SortedByCoordinate
stringtie SRR500121Aligned.sortedByCoord.out.bam -A SRR500121.tab -o SRR500121.gtf -G gencode.v24.annotation.gtf -m 100 -c 0.001 -e


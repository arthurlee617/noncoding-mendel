conda activate tobias
 
TOBIAS ATACorrect --bam /home/tobias/CN7-10-5/CN7-e10-5_merged.bam --genome /home/arthur/tobias/GRCm38.p4.genome.fa --peaks /home/tobias/CN7-10-5/CN7_e10-5_mergedPeaks.bed --blacklist /home/tobias/mm10-blacklist.v2.bed --outdir --cores 60

TOBIAS ATACorrect --bam /home/tobias/CN7-11-5/CN7-e11-5_merged.bam --genome /home/tobias/GRCm38.p4.genome.fa --peaks /home/tobias/CN7-11-5/CN7_e11-5_mergedPeaks.bed --blacklist /home/tobias/mm10-blacklist.v2.bed --outdir /home/tobias/CN7-11-5 --cores 60

TOBIAS FootprintScores --signal /home/tobias/corrected/CN7-e10-5_merged_corrected.bw --regions /home/tobias/CN7-10-5/CN7_e10-5_mergedPeaks.bed --output CN7_10-5_footprints.bw --cores 60

TOBIAS ScoreBed --bed /home/tobias/CN7-10-5/CN7_e10-5_mergedPeaks.bed --bigwigs /home/tobias/corrected/CN7-e10-5_merged_corrected.bw --output /home/tobias/CN7-10-5/merged_peaks_scored.bed

TOBIAS TFBScan --motifs /home/tobias/motifs.jaspar --fasta /home/arthur/tobias/GRCm38.p4.genome.fa --cores 60 --regions /home/tobias/CN7-10-5/CN7_e10-5_mergedPeaks.bed --add-region-columns --outdir /home/tobias/CN7-10-5/tfbscan_output

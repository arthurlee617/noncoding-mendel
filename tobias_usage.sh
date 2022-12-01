conda activate tobias
 
TOBIAS ATACorrect --bam /home/arthur/tobias/CN7-10-5/CN7-e10-5_merged.bam --genome /home/arthur/tobias/GRCm38.p4.genome.fa --peaks /home/arthur/tobias/CN7-10-5/CN7_e10-5_mergedPeaks.bed --blacklist /home/arthur/tobias/mm10-blacklist.v2.bed --outdir corrected --cores 60

TOBIAS ATACorrect --bam /home/arthur/tobias/CN7-11-5/CN7-e11-5_merged.bam --genome /home/arthur/tobias/GRCm38.p4.genome.fa --peaks /home/arthur/tobias/CN7-11-5/CN7_e11-5_mergedPeaks.bed --blacklist /home/arthur/tobias/mm10-blacklist.v2.bed --outdir /home/arthur/tobias/CN7-11-5 --cores 60

TOBIAS FootprintScores --signal /home/arthur/tobias/corrected/CN7-e10-5_merged_corrected.bw --regions /home/arthur/tobias/CN7-10-5/CN7_e10-5_mergedPeaks.bed --output CN7_10-5_footprints.bw --cores 60

TOBIAS ScoreBed --bed /home/arthur/tobias/CN7-10-5/CN7_e10-5_mergedPeaks.bed --bigwigs /home/arthur/tobias/corrected/CN7-e10-5_merged_corrected.bw --output /home/arthur/tobias/CN7-10-5/merged_peaks_scored.bed

TOBIAS TFBScan --motifs /home/arthur/tobias/motifs.jaspar --fasta /home/arthur/tobias/GRCm38.p4.genome.fa --cores 60 --regions /home/arthur/tobias/CN7-10-5/CN7_e10-5_mergedPeaks.bed --add-region-columns --outdir /home/arthur/tobias/CN7-10-5/tfbscan_output

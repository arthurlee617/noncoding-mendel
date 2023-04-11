conda activate basenji

/home/arthur/basenji/bin/basenji_sad.py -f /home/arthur/lauren/basenji_training/data/fasta/GRCh38.primary_assembly.genome.fa -o /home/arthur/sandbox/sadness/chr3_sad_2023 --rc --shift "1,0,-1" -t /home/arthur/lauren/basenji_training/data/cn-wigs.txt /home/arthur/lauren/basenji_training/models/cn/params.json /home/arthur/lauren/basenji_training/models/cn/model_best.h5 /home/arthur/sandbox/sadness/chr3_2023.vcf 

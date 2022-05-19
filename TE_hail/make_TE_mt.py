import hail as hl
hl.init()

mt = hl.import_matrix_table('TE-900Genomes.tsv', no_header = False, row_fields={'Chr':hl.tstr,'Gene_ann':hl.tstr,'X1KGP_ann':hl.tstr,'X1KGP_AF':hl.tstr,'X1KGP_strand':hl.tstr,'EuL1db_ann':hl.tstr,'EuL1db_AF':hl.tstr,'EuL1db_strand':hl.tstr,'Evrony_ann':hl.tstr,'Evrony_strand':hl.tstr,'KNR_ann':hl.tstr,'KNR_papers':hl.tstr,'Overlapped_repeats':hl.tstr,'Overlapped_homopolymers':hl.tstr,'gnomAD_ann':hl.tstr,'gnomAD_AF':hl.tstr,'gnomAD_strand':hl.tstr,'gnomAD_pLI':hl.tstr,'locus':hl.tstr, 'TE.type':hl.tstr}, entry_type=hl.tstr, row_key = 'locus') 
#Import the TE calls as a matrix table, ref genome hg19
#Fix formatting of sample ID's to match 900-genomes; - vs .

mt.aggregate_rows(hl.agg.collect_as_set(mt.Chr)) #check for bad chromosomes

bad_chroms = ["Un_gl000220","1_gl000192_random","Un_gl000212","4_gl000193_random","Un_gl000241", "Un_gl000237", '17_gl000203_random','18_gl000207_random','19_gl000208_random','1_gl000192_random','4_gl000193_random','8_gl000196_random','8_gl000196_random','9_gl000198_random',
 '9_gl000199_random','Un_gl000211','Un_gl000212','Un_gl000217','Un_gl000220','Un_gl000229', 'Un_gl000231','Un_gl000232','Un_gl000233','Un_gl000234',
 'Un_gl000235','Un_gl000237','Un_gl000241'] 
bad_chroms = hl.literal(bad_chroms)
mt = mt.filter_rows(bad_chroms.contains(mt.Chr), keep=False) #filter out bad chromosomes

mt = mt.annotate_rows(new_key = hl.parse_locus(mt.locus, reference_genome='GRCh37')) #Convert to locus from string
mt = mt.key_rows_by(mt.new_key)
mt = mt.drop('locus')
mt = mt.annotate_rows(locus = mt.new_key)
mt = mt.key_rows_by(mt.locus)
mt = mt.drop('new_key')
mt = mt.annotate_cols(s = mt.col_id)
mt = mt.key_cols_by('s')
mt = mt.drop('col_id')

rg37 = hl.get_reference('GRCh37')
rg38 = hl.get_reference('GRCh38')  
rg37.add_liftover('gs://hail-common/references/grch37_to_grch38.over.chain.gz', rg38) #from gs://hail-common/references/grch37_to_grch38.over.chain.gz 

mt = mt.annotate_rows(new_locus=hl.liftover(mt.locus, 'GRCh38', include_strand=True),
                  old_locus=mt.locus)  #liftover, keep old locus 
mt = mt.filter_rows(hl.is_defined(mt.new_locus) & ~mt.new_locus.is_negative_strand)  
mt = mt.key_rows_by(locus=mt.new_locus.result) 

mt.write('TE-all.mt', overwrite = True) 
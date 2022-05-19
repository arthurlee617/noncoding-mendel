import hail as hl
hl.init()
TE_mt = hl.read_matrix_table('TE-all.mt')
TE_mt = TE_mt.drop('Chr',
    'Gene_ann',
    'X1KGP_ann',
    'X1KGP_AF',
    'X1KGP_strand',
    'EuL1db_ann',
    'EuL1db_AF',
    'EuL1db_strand',
    'Evrony_ann',
    'Evrony_strand',
    'KNR_ann',
    'KNR_papers',
    'Overlapped_repeats',
    'Overlapped_homopolymers',
    'gnomAD_ann',
    'gnomAD_AF',
    'gnomAD_strand',
    'gnomAD_pLI',
    'new_locus',
    'old_locus')

TE_mt = TE_mt.annotate_rows(alleles = hl.empty_array(hl.tstr))
TE_mt = TE_mt.annotate_entries(GT = hl.null('call'))
TE_mt = TE_mt.drop('TE.type')

mt = hl.read_matrix_table('G900-filteredregions-variants.mt')

mt_SV = mt.key_rows_by('locus')
mt_SV = mt_SV.drop('rsid', 'qual', 'filters', 'info', 'AD','DP', 'GQ', 'MIN_DP','PGT','PID','PL','RGQ','SB','AB','a_index','was_split','topmed_freq','gnomad_genomes_freq','gerp_scores', 'gnomad_exomes_freq',
 'gnomad_genomes_AC',
 'gnomad_exomes_AC',
 'variant_qc',
 'most_severe_consequence',
 'clinvar_variant_summary',
  'peak_label', 'peds')



mt_SV = mt_SV.annotate_entries(x = hl.null('str'))

mt_SV = mt_SV.key_cols_by('s')
TE_mt = TE_mt.key_cols_by('s')

mt_SV = mt_SV.rename({'GT': 'GT2', 'x': 'x2'})
mt_SV = mt_SV.select_entries(x=mt_SV.x2, GT=mt_SV.GT2)

union_MT = mt_SV.union_rows(TE_mt)

union_MT.write("TE-shortvar-union.mt", overwrite=True)
import hail as hl
import argparse
from gnomad.utils.vep import vep_or_lookup_vep, process_consequences
from gnomad.resources.grch38 import gnomad
YOUR_GCP_PROJECT_NAME = 'cmg-test'
hl.init(default_reference='GRCh38',spark_conf={'spark.hadoop.fs.gs.requester.pays.mode': 'AUTO',
    'spark.hadoop.fs.gs.requester.pays.project.id': YOUR_GCP_PROJECT_NAME})
gnomad_v3_genomes = gnomad.public_release("genomes")
gnomad_genomes = gnomad_v3_genomes.ht()
gnomad_exomes = hl.read_table('gs://gnomad-public-requester-pays/release/2.1.1/liftover_grch38/ht/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.ht') 
rg = hl.get_reference('GRCh38')
#gnomad_genomes.write('gnomad.ht')
def parse_args():
    p = argparse.ArgumentParser(description="Create and write to disk annotation table from hail annotation DB on gcp for running hail queries locally")

    p.add_argument("input_mt_path", help="path of input matrix table (.mt) containing variants to annotate")
    p.add_argument("-f", "--fam-file", help="path of .fam file that encodes pedigree relationships between samples in "
                   "the dataset. --family-id and --inheritance args can't be used unless --fam-file is specified.")
    
    p.add_argument("-o", "--output-ht", help="path of output .ht file.")
    args = p.parse_args()

    return args
def custom_split_multi_hts(ds, keep_star=False, left_aligned=False, vep_root='vep', *, permit_shuffle=False):
    '''This function adapted from hail source code for hl.split_multi_hts to split our 900 genomes VCF, with vcf imported with PGT/PID tarray field types.
    Splits multi-allelic variants. '''
    split = hl.split_multi(ds, keep_star=keep_star, left_aligned=left_aligned, permit_shuffle=permit_shuffle)
    row_fields = set(ds.row)
    update_rows_expression = {}
    if vep_root in row_fields:
        update_rows_expression[vep_root] = split[vep_root].annotate(**{
            x: split[vep_root][x].filter(lambda csq: csq.allele_num == split.a_index)
            for x in ('intergenic_consequences', 'motif_feature_consequences',
                      'regulatory_feature_consequences', 'transcript_consequences')})

    if isinstance(ds, hl.Table):
        return split.annotate(**update_rows_expression).drop('old_locus', 'old_alleles')

    split = split.annotate_rows(**update_rows_expression)
    entry_fields = ds.entry

    expected_field_types = {
        'GT': hl.tcall,
        'AD': hl.tarray(hl.tint),
        'DP': hl.tint,
        'GQ': hl.tint,
        'PL': hl.tarray(hl.tint),
        'PGT': hl.tarray(hl.tstr),
        'PID': hl.tarray(hl.tstr)
    }

    bad_fields = []
    for field in entry_fields:
        if field in expected_field_types and entry_fields[field].dtype != expected_field_types[field]:
            bad_fields.append((field, entry_fields[field].dtype, expected_field_types[field]))

    if bad_fields:
        msg = '\n  '.join([f"'{x[0]}'\tfound: {x[1]}\texpected: {x[2]}" for x in bad_fields])
        raise TypeError("'split_multi_hts': Found invalid types for the following fields:\n  " + msg)

    update_entries_expression = {}
    if 'GT' in entry_fields:
        update_entries_expression['GT'] = hl.downcode(split.GT, split.a_index)
    if 'DP' in entry_fields:
        update_entries_expression['DP'] = split.DP
    if 'AD' in entry_fields:
        update_entries_expression['AD'] = hl.or_missing(hl.is_defined(split.AD),
                                                        [hl.sum(split.AD) - split.AD[split.a_index], split.AD[split.a_index]])
    if 'PL' in entry_fields:
        pl = hl.or_missing(
            hl.is_defined(split.PL),
            (hl.range(0, 3).map(lambda i:
                                hl.min((hl.range(0, hl.triangle(split.old_alleles.length()))
                                        .filter(lambda j: hl.downcode(hl.unphased_diploid_gt_index_call(j),
                                                                      split.a_index).unphased_diploid_gt_index() == i
                                                ).map(lambda j: split.PL[j]))))))
        if 'GQ' in entry_fields:
            update_entries_expression['PL'] = pl
            update_entries_expression['GQ'] = hl.or_else(hl.gq_from_pl(pl), split.GQ)
        else:
            update_entries_expression['PL'] = pl
    else:
        if 'GQ' in entry_fields:
            update_entries_expression['GQ'] = split.GQ

    if 'PGT' in entry_fields:
        update_entries_expression['PGT'] = split.PGT
    if 'PID' in entry_fields:
        update_entries_expression['PID'] = split.PID
    return split.annotate_entries(**update_entries_expression).drop('old_locus', 'old_alleles')

def main():
    args = parse_args()
    mt = hl.read_matrix_table(args.input_mt_path)
    mt = hl.split_multi_hts(mt, keep_star = False) #split multiallelic variants
    # Annotate with VEP: most severe consequence
   # ds = vep_or_lookup_vep(mt.rows(), reference="GRCh38")
    #ds = process_consequences(ds)
    #mt = mt.annotate_rows(most_severe_consequence = ds[mt.locus, mt.alleles].vep.most_severe_consequence)
    #Annotate with pedigree info
    if args.fam_file:
        ped = hl.import_table(args.fam_file, no_header=True)
        ped = ped.select(famID = ped.f0, proband_ID = ped.f1, father_ID = ped.f2, mother_id= ped.f3, sex = ped.f4, affected=ped.f5)
        ped = ped.key_by('proband_ID')
        mt = mt.annotate_cols(peds=ped[mt.s]) # annotate with pedigree 
    mt = mt.annotate_entries(AB = hl.min(mt.AD[0]/hl.sum(mt.AD), mt.AD[1]/hl.sum(mt.AD))) # AD: [Ref, Alt1, Alt2, etc]
    mt = mt.annotate_rows(gnomad_genomes_freq = hl.if_else(hl.is_defined(gnomad_genomes[mt.locus, mt.alleles].freq.AF[1]),gnomad_genomes[mt.locus, mt.alleles].freq.AF[1],0)) #annotate with gnomad genomes freq
    mt = mt.annotate_rows(gnomad_genomes_AC = hl.if_else(hl.is_defined(gnomad_genomes[mt.locus, mt.alleles].freq.AC[1]),gnomad_genomes[mt.locus, mt.alleles].freq.AC[1],0)) 
    mt = mt.annotate_rows(gnomad_exomes_freq = hl.if_else(hl.is_defined(gnomad_exomes[mt.locus, mt.alleles].freq.AF[1]),gnomad_exomes[mt.locus, mt.alleles].freq.AF[1],0)) #annotate with gnomad exomes freq
    mt = mt.annotate_rows(gnomad_exomes_AC = hl.if_else(hl.is_defined(gnomad_exomes[mt.locus, mt.alleles].freq.AC[1]),gnomad_exomes[mt.locus, mt.alleles].freq.AC[1],0)) 
    mt = hl.variant_qc(mt)
   # gerp = hl.read_table('gs://hail-datasets-us/annotations/GERP_scores.GERP++.GRCh38.ht/')
    #gerp.write('gerpscores.ht')
    #mt = hl.vep(mt)
   # mt = db.annotate_rows_db(mt,  "gerp_scores", "gerp_elements") #annotate with clinvar
    mt.rows().write(args.output_ht)

if __name__ == "__main__":
    main()


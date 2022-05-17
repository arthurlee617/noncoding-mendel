import argparse
import hail as hl
import logging
import os
import re

hl.init(default_reference='GRCh38')


logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s')
logger = logging.getLogger()
logger.setLevel(logging.INFO)

def parse_args():
    p = argparse.ArgumentParser(description="This script lets you filter variants from a hail .mt (matrix table) based "
        "on their annotations and/or genotypes and outputs the passing variants to a .tsv (tab-separated) text file. "
        "All arguments are optional except the path of the input matrix table. If no other args are specified, all "
        "variants will be written to the output .tsv file up to the limit specified by -M (which defaults to 100,000 "
        "variants). If one or more filter args is specified, each variant will only be included in the output file if "
        "it passes all the specified filters. In other words, the filters are applied serially. On the other hand, "
        "any filter args that need a non-numerical value (eg. --gene-id) also support multiple values separated by "
        "commas. The commas are interpreted as logical OR's. For example, --gene-id ENSG00000198947,ENSG00000183091 "
        "will output variants that are in ENSG00000198947 OR in ENSG00000183091.")

    p.add_argument("-f", "--fam-file", help="path of .fam file that encodes pedigree relationships between samples in "
                   "the dataset. --family-id and --inheritance args can't be used unless --fam-file is specified.")


    p.add_argument("-i", "--family-id", action="append", help="a family id found in the 1st column of the .fam "
                   "file. This arg and --fam-file must be specified before --inheritance filter can be used. "
                   "--family-id can be specified more than once, in which case variants will only pass "
                   "an --inheritance filter if they pass it in all the specified families.")

    p.add_argument("-H", "--inheritance", help="filter to variants whose genotypes are consistent with a particular "
                   "inheritance mode in the family(s) specified by the --family-id arg. Multiple inheritance values "
                   "can be provided, separated by commas (,) meaning 'OR'. Comp-het is not active yet. "")")

    p.add_argument("-s", "--strict", action="store_true", help="modifies the behavior of the "
                   "--inheritance filter to be more strict. Without this, inheritance filters allow for inaccurate "
                   "genotype calls. By default, the 'dominant' filter allows both parents to be homozygous reference "
                   "in case one of their genotypes is actually het. In that case, specifying -s enforces that exactly "
                   "one parent should be het. Similarly, by default, recessive filters allow one or both parents to be "
                   "homozygous reference, and specifying -s enforces that both parents should be het.") ## not functional ? 


    p.add_argument("-m", "--allow-missing-genotypes", action="store_true", help="modifies the behavior of the "
                   "--inheritance filter to allow missing genotypes. -m and -s can be specified together, in which "
                   "case --inheritance filters will allow missing genotypes, but not homozygous reference genotypes "
                   "(see description for -s).")

    p.add_argument("-gp", "--gerpscore", type=float, help="Minimum GERP score.")

    p.add_argument("--annotations", help="filter by VEP transcript consequence (https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html), clinvar pathogenicity or HGMD disease. "
                   "Multiple values can be provided, separated by commas (,) or pluses (+) meaning 'OR'.", choices=[
                        "clinvar-p", "clinvar-lp", "clinvar-vus",
                        
                           "frameshift_variant", "missense_variant", "coding_sequence_variant","TF_binding_site_variant","regulatory_region_amplification","regulatory_region_variant"
                    ]) 

    p.add_argument("--gnomad-genomes-af", type=float, help="Maximum global allele frequency in gnomAD genomes")
    p.add_argument("--gnomad-genomes-ac", type=float, help="Maximum global allele count in gnomAD genomes")
    #GERP
    #Conservation
    p.add_argument("--gnomad-exomes-af", type=float, help="Maximum global allele frequency in gnomAD exomes")
    p.add_argument("--gnomad-exomes-ac", type=float, help="Maximum global allele count in gnomAD exomes")
    p.add_argument("--topmed-af", type=float, help="Maximum global allele frequency in TopMed")
    p.add_argument("--topmed-ac", type=float, help="Maximum global allele count in TopMed")
    p.add_argument("--this-callset-af", type=float, help="Maximum allele frequency in the input callset")
    p.add_argument("--this-callset-ac", type=float, help="Maximum allele count in the input callset")

    p.add_argument("-g", "--genes", help="ENSG ensembl id or gene name. Only return variants in the given gene(s). "
                   "Multiple values can be provided, separated by commas. Alternatively, the value can be the path of "
                   "a text file which contains one or more genes.") #annotated with VEP 

    p.add_argument("-L", "--intervals", help="<chrom>:<start>-<end> Only return variants in this genomic region. "
                   "Multiple values can be provided, separated by commas. Alternatively, the value can be the path of "
                   "a text file which contains one or more such intervals.")

    p.add_argument("--pass-filter", action="store_true", help="only output variants that have 'PASS' in the VCF "
                   "filter column")
    p.add_argument("--filter-unaffected", action="store_true", help="only output variants that have high affected ratio"
                   )
    p.add_argument("--filter-rec", action="store_true", help="only output variants that do not have unaffected 1/1 individuals"
                   )

    p.add_argument("-gq", "--gq", type=int, help="minimum genotype quality (GQ) for all samples in families "
                   "specified by --family-id")
    p.add_argument("-ab", "--ab", type=float, help="minimum allele balance (AB) for all samples in families "
                   "specified by --family-id. This should be a decimal between 0 and 1")
    p.add_argument("--hail-db", help="Path to local annotations table"
                   )

    p.add_argument("-gt", "--gt", action="append", help="sample genotypes filter. Only output variants where a "
                   "specific sample has the given genotype(s). The value have the format <sample-id>=<genotype>"
                   " where <genotype> can be HOM-REF, HET, HOM-ALT, MISSING, or a comma-separated list of these. "
                   "For example, --gt F32-M=HET,MISSING means only output variants where the genotype of "
                   "sample 'F32-M' was either HET (0/1) or MISSING (./.). This argument can be specified "
                   "more than once to apply a genotype filter to more than one sample id. For example, "
                   "--gt F32-M=HET,MISSING --gt F32-P=HOM would only output variants that pass both genotype filters.")

    p.add_argument("-w", "--write-file",  help="Write filtered MT to file to reduce computation time"
    				"(Improves run time if interval filter used)")

    p.add_argument("-o", "--output-tsv", help="path of output .tsv file.")
    p.add_argument("-M", "--max-variants", help="the maximum number of variants to write to the output .tsv file. "
                   "This is an optional circuit breaker for preventing the output file from being too large.",
                   default=100_000)

    p.add_argument("input_mt_path", help="path of input matrix table (.mt)")


    args = p.parse_args()

    if args.strict: p.error("--strict arg isn't implemented yet")
    if args.allow_missing_genotypes: p.error("--allow-missing-genotypes arg isn't implemented yet")

    return args


INTERVAL_REGEXP = re.compile("((chr)?[0-9XYMT]{1,2})[:\s]+([0-9]{1,9})[-\s]+([0-9]{1,9})", re.IGNORECASE)

### pass-filter
#rg = hl.get_reference('GRCh38')
#rg.add_sequence('gs://hail-common/references/Homo_sapiens_assembly38.fasta.gz','gs://hail-common/references/Homo_sapiens_assembly38.fasta.fai')

def custom_split_multi_hts(ds, keep_star=False, left_aligned=False, vep_root='vep', *, permit_shuffle=False):
    '''adapted from hail source code for hl.split_multi_hts to split our MT with vcf imported with PGT/PID tarray field types.
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

def filter_by_ensembl_geneID(mt, genes):
    # Parsing VEP annotations is very tricky! This works. 
    mt = mt.annotate_rows(gene_id = mt.vep.transcript_consequences.gene_ID)
    genes_table = mt.gene_id 
    filtered_genes = genes_table.map(lambda x: hl.array([x.contains(gene) for gene in genes]).contains(True))
    mt = mt.annotate_rows(filtered = filtered_genes) #annotate mt with this table
    mt = mt.filter_rows(hl.is_defined(mt.filtered), keep = True) #keep only rows with t/f annotations
    mt = mt.filter_rows(mt.filtered.contains(True)) # keep only rows with t annotations 
    return mt
def filter_by_geneSymbol(mt, genes): 
    mt = mt.annotate_rows(gene_id = mt.vep.transcript_consequences.gene_symbol)
    genes_table = mt.gene_id 
    filtered_genes = genes_table.map(lambda x: hl.array([x.contains(gene) for gene in genes]).contains(True))
    mt = mt.annotate_rows(filtered = filtered_genes) #annotate mt with this table
    mt = mt.filter_rows(hl.is_defined(mt.filtered), keep = True) #keep only rows with t/f annotations
    mt = mt.filter_rows(mt.filtered.contains(True)) # keep only rows with t annotations 
    return mt
def filter_by_genes(mt, genes):
    #this assumes file is separated by , tab, newline or whitespace
    if os.path.isfile(genes):
        with open(genes, "rt") as f:
            genes = f.read()
            genes_list = re.split(',|\t|\n| ', genes)
    elif genes.startswith("gs://"):
        with open(genes, "rt") as f:
            genes = f.read()
            genes_list = re.split(',|\t|\n| ', genes)
    else:
        genes_list = genes.split(',')
    geneID = []
    genesymbol = []
    for gene in genes_list:
        if gene.startswith('ENSG'):
            geneID.append(gene)
        else:
            genesymbol.append(gene)
    if geneID:
        mt = filter_by_ensembl_geneID(mt, geneID)
    if genesymbol:
        mt = filter_by_geneSymbol(mt, genesymbol)
    return mt

def apply_topmed_AF(mt, topmedAF):
    mt = mt.filter_rows((mt.topmed_freq  < topmedAF) | hl.is_missing(mt.topmed_freq), keep = True) # [1] index based on: https://broadinstitute.github.io/gnomad_methods/_modules/gnomad/utils/filtering.html#filter_by_frequency
    # Descriptions for this list in gnomad_genomes.freq_meta [1] -> {'group': 'raw'} (0: {'group': 'adj'}; latter are stratified populations)
    return mt

def apply_max_gnomad_genomes_AF(mt, ggafARG):
    mt = mt.filter_rows(mt.gnomad_genomes_freq  < ggafARG , keep = True) # [1] index based on: https://broadinstitute.github.io/gnomad_methods/_modules/gnomad/utils/filtering.html#filter_by_frequency
    # Descriptions for this list in gnomad_genomes.freq_meta [1] -> {'group': 'raw'} (0: {'group': 'adj'}; latter are stratified populations)
    return mt
def apply_max_gnomad_exomes_AF(mt, geafARG):
    mt = mt.filter_rows(mt.gnomad_exomes_freq < geafARG, keep = True)
    return mt
def apply_max_gnomad_genomes_AC(mt, ggacARG):  
    mt = mt.filter_rows(mt.gnomad_genomes_AC < ggacARG , keep = True) # [1] index based on: https://broadinstitute.github.io/gnomad_methods/_modules/gnomad/utils/filtering.html#filter_by_frequency
    # Descriptions for this list in gnomad_genomes.freq_meta [1] -> {'group': 'raw'} (0: {'group': 'adj'}; latter are stratified populations)
    return mt
def apply_max_gnomad_exomes_AC(mt, geacARG):
    mt = mt.filter_rows(mt.gnomad_exomes_AC  < geacARG , keep = True)
    return mt
def pass_filter(mt):
    mt = mt.filter_rows(mt.filters.size() > 0, keep=False)
    return mt
# Annotations - Clinvar ONLY so far (8/12)
def apply_annotations_filter(mt, annotation):
    annotations = annotation.split(",")
    mt = mt.annotate_rows(clinsig_unique_annotations = hl.set(mt.clinvar_variant_summary['ClinicalSignificance']).length())
    mt = mt.annotate_rows(UnanimousClinVar = mt.clinvar_variant_summary['ClinicalSignificance'][0])
    conflicting = mt.clinsig_unique_annotations > 1 # conflicting conditions
    vus = (mt.UnanimousClinVar =='Uncertain significance')
    p = ((mt.UnanimousClinVar =='Pathogenic') & (mt.clinsig_unique_annotations == 1))
    lp = ((mt.UnanimousClinVar =='Likely pathogenic')& (mt.clinsig_unique_annotations == 1))
    lb = ((mt.UnanimousClinVar =='Likely benign')& (mt.clinsig_unique_annotations == 1))
    b = ((mt.UnanimousClinVar == 'Benign') & (mt.clinsig_unique_annotations == 1))
    conds = False
    if "clinvar-vus" in annotations:
        conds = conds | vus | conflicting
    if "clinvar-p" in annotations: 
        conds = conds | p
    if "clinvar-lp" in annotations: 
        conds = conds | lp
    if "clinvar-lb" in annotations:
        conds = conds | lb
    if "clinvar-b" in annotations: 
        conds = conds | b
    transcript_cons = False
    for a in annotations:
        transcript_cons = transcript_cons | (mt.most_severe_consequence == a)

    mt = mt.filter_rows((hl.agg.count_where(conds) >  0) | transcript_cons, keep = True)

    return mt

#def apply_transcript_consequences(mt, consequences):
    ##mt.vep.most_severe_consequence
    ##https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html
### GT
def sample_genotype_conditions(mt, sample, gt):
    if gt == 'HOM-REF': #return MT with variants where sample is hom-ref
        condition = hl.agg.count_where((mt.s == sample) & mt.GT.is_hom_ref()) > 0
    elif gt == 'HET': #return MT with varints where sample is het
        condition = hl.agg.count_where((mt.s == sample) & mt.GT.is_het()) > 0
    elif gt == 'HOM-ALT': #return MT with variants where sample is hom-var
        condition = hl.agg.count_where((mt.s == sample) & mt.GT.is_hom_var()) > 0 
    elif gt == 'MISSING': #Return MT with variants where sample GT is missing
        condition = hl.agg.count_where((mt.s == sample) & (hl.is_defined(mt.GT) == False)) > 0
    else:
        print('Not valid GT')
    mt = mt.filter_rows(condition, keep=True)
    return mt
def apply_sample_genotype(mt, gtarg):
    #['sample=gt,gt,gt', 'sample=gt,gt,gt']
    # Initialize return MT by getting first sample-GT combo, creat MT based on filtered conditions
    sample, gt = (gtarg[0].split('=')[0], gtarg[0].split('=')[1].split(',')[0])
    #mt = sample_genotype_conditions(mt, sample, gt)
    #mt = mt.key_rows_by(mt['locus'], mt['alleles'])
    for arg in gtarg:
        # iterate through all other combos
        sample = arg.split('=')[0]
        genotypes = arg.split('=')[1].split(',')
        for gt in genotypes:
            mt = sample_genotype_conditions(mt, sample, gt)
            mt = mt.key_rows_by(mt['locus'], mt['alleles'])
    return mt
### GQ
def apply_min_GQ(mt, GQarg):
   ## this computes the min per row. see https://hail.is/docs/0.2/hail.MatrixTable.html#hail.MatrixTable.annotate_rows
    mt = mt.filter_entries(mt.GQ > GQarg, keep = True)
    return mt
### Allele balance
# AD stored as [X,X]
# Need to split multiallelic variants before doing this
def apply_AB_filter(mt, minABarg):
    mt = mt.filter_entries(mt.GT.is_het() & (mt.AB < minABarg), keep=False) # filter out all entries that are het and have an allele balance of less than the threshold
    return mt

def apply_inher_filter(mt, famFile, famIDarg, inherArg):
    inherModes = inherArg.split(",")
    trimmed_mt = mt
    pedigree = hl.Pedigree.read(famFile)
    deNovo = hl.de_novo(mt, pedigree, hl.cond(hl.is_defined(mt.gnomad_genomes_freq), mt.gnomad_genomes_freq, mt.variant_qc.AF[1]))
    rows_to_keep = family_pass_vars(trimmed_mt, famFile, famIDarg, inherModes, deNovo).rows()
    rows_to_keep = rows_to_keep.key_by('locus','alleles')
    mt = mt.semi_join_rows(rows_to_keep)
    return mt

def family_pass_vars(mt, famFile, famIDarg, inherModes, deNovo):
    #Takes MT annotated with pedigree, fam file, singular family ID, singular inheritance mode
    pedigree = hl.Pedigree.read(famFile)
    trio_mt = hl.trio_matrix(mt, pedigree, complete_trios=True) #Until we can deal with missing genotypes, keep complete trios only
    to_keep = hl.literal(famIDarg)
    trio_mt = trio_mt.filter_cols(to_keep.contains(trio_mt.fam_id))
    trio_mt = trio_mt.filter_cols(trio_mt.proband['peds']['affected'] =='2')# Trio inheritance cases where proband is affected
    mother_affected = trio_mt.mother['peds']['affected'] =='2'
    father_affected = trio_mt.father['peds']['affected'] =='2'
    mother_unaffected = trio_mt.mother['peds']['affected'] =='1'
    father_unaffected = trio_mt.father['peds']['affected'] =='1'


    MorFhet_dom = (((trio_mt.father_entry['GT'].is_het() & father_affected & trio_mt.mother_entry['GT'].is_hom_ref() & mother_unaffected ) 
                | (trio_mt.mother_entry['GT'].is_het() & trio_mt.father_entry['GT'].is_hom_ref() & father_unaffected & mother_affected)) 
               & trio_mt.proband_entry['GT'].is_het()) #One het/one hom ref parent. Het parent affected, child is het and affected. 


    MandFhet_hetProband_dom = ((trio_mt.father_entry['GT'].is_het() & father_affected & trio_mt.mother_entry['GT'].is_het() & mother_affected) 
               & trio_mt.proband_entry['GT'].is_het() ) #Both parents are het and affected, child is het and affected

    MandFhet_homProband_rec = ((trio_mt.father_entry['GT'].is_het() & father_unaffected & trio_mt.mother_entry['GT'].is_het() & mother_unaffected) 
               & trio_mt.proband_entry['GT'].is_hom_var() ) #Both parents are het and unaffected , child is hom-non-ref and affected

    MandFhet_homProband_dom = ((trio_mt.father_entry['GT'].is_het() & father_affected & trio_mt.mother_entry['GT'].is_het() & mother_affected) 
               & trio_mt.proband_entry['GT'].is_hom_var() ) #Both parents are het affected , child is hom-non-ref and affected

    MorFhet_homNR_dom = (((trio_mt.father_entry['GT'].is_het() & father_affected & trio_mt.mother_entry['GT'].is_hom_var() & mother_affected) 
                | (trio_mt.mother_entry['GT'].is_het() & mother_affected & trio_mt.father_entry['GT'].is_hom_var() & father_affected)) 
               & trio_mt.proband_entry['GT'].is_het()) #One het/one hom var, both parents affected, child is het and affected

    MorFhet_homNR_homProband_dom = (((trio_mt.father_entry['GT'].is_het() & father_affected & trio_mt.mother_entry['GT'].is_hom_var() & mother_affected) 
                | (trio_mt.mother_entry['GT'].is_het() & father_affected & mother_affected & trio_mt.father_entry['GT'].is_hom_var())) 
               & trio_mt.proband_entry['GT'].is_hom_var()) #One het/one hom var, child is hom, all affected

    MorFhet_homNR_homProband_rec = (((trio_mt.father_entry['GT'].is_het() &father_unaffected & trio_mt.mother_entry['GT'].is_hom_var() &mother_affected) 
                | (trio_mt.mother_entry['GT'].is_het() & mother_unaffected & trio_mt.father_entry['GT'].is_hom_var() &father_affected)) 
               & trio_mt.proband_entry['GT'].is_hom_var()) #One het/one hom var (unaffected/affected, child is hom and affected


    MandFhomvar_domrec = (((trio_mt.father_entry['GT'].is_hom_var() & father_affected& trio_mt.mother_entry['GT'].is_hom_var() &mother_affected)) 
               & trio_mt.proband_entry['GT'].is_hom_var()) #all three are 1/1, all affected (dom or recessive)

    MorFhomvar_dom = (((trio_mt.father_entry['GT'].is_hom_var() & father_affected & trio_mt.mother_entry['GT'].is_hom_ref() &mother_unaffected) |
                    (trio_mt.father_entry['GT'].is_hom_ref() & father_unaffected& trio_mt.mother_entry['GT'].is_hom_var() &mother_affected) ) 
               & trio_mt.proband_entry['GT'].is_het()) #one parent 1/1, one parent 0/0, child 0/1, dominant

    MorFhet_rec = (((trio_mt.father_entry['GT'].is_het() & father_unaffected & trio_mt.mother_entry['GT'].is_hom_ref() &mother_unaffected) |
                    (trio_mt.father_entry['GT'].is_hom_ref() & father_unaffected& trio_mt.mother_entry['GT'].is_het() &mother_unaffected) ) 
               & trio_mt.proband_entry['GT'].is_hom_var()) #one parent 0/1, one parent 0/0, child 1/1, rec 

    ## Xlinked
    Xlinked1 = (((trio_mt.mother_entry['GT'].is_hom_ref() & trio_mt.father_entry['GT'].is_non_ref() & mother_unaffected & father_affected)) 
               & trio_mt.proband_entry['GT'].is_non_ref() & (trio_mt.proband.peds.sex == '1'))

    Xlinked2 = (((trio_mt.mother_entry['GT'].is_het() & trio_mt.father_entry['GT'].is_hom_ref() & mother_unaffected & father_unaffected)) 
               & trio_mt.proband_entry['GT'].is_non_ref() & (trio_mt.proband.peds.sex == '1'))

    Xlinked3 = (((trio_mt.mother_entry['GT'].is_het() & trio_mt.father_entry['GT'].is_non_ref() & mother_unaffected & father_affected)) 
               & trio_mt.proband_entry['GT'].is_hom_var() & (trio_mt.proband.peds.sex == '2'))

    Xlinked4 = (((trio_mt.mother_entry['GT'].is_het() & trio_mt.father_entry['GT'].is_non_ref() & mother_unaffected & father_affected)) 
               & trio_mt.proband_entry['GT'].is_non_ref() & (trio_mt.proband.peds.sex == '1'))

    Xlinked5 = (((trio_mt.mother_entry['GT'].is_hom_var() & trio_mt.father_entry['GT'].is_non_ref() & mother_affected & father_affected)) 
               & trio_mt.proband_entry['GT'].is_non_ref() & (trio_mt.proband.peds.sex == '1'))

    Xlinked6 = (((trio_mt.mother_entry['GT'].is_hom_var() & trio_mt.father_entry['GT'].is_non_ref() & mother_affected & father_affected)) 
               & trio_mt.proband_entry['GT'].is_hom_var() & (trio_mt.proband.peds.sex == '2'))
    
    X_linked_conditions =  hl.agg.count_where((Xlinked1 | Xlinked2 | Xlinked3 | Xlinked4 | Xlinked5 | Xlinked6 )) > 0
    dominant = hl.agg.count_where((MorFhet_dom | MandFhet_hetProband_dom | MandFhet_homProband_dom | MorFhet_homNR_dom | MorFhet_homNR_homProband_dom | MandFhomvar_domrec | MorFhomvar_dom)) > 0
    homrecessive = hl.agg.count_where((MandFhet_homProband_rec | MorFhet_homNR_homProband_rec | MandFhomvar_domrec | MorFhet_rec)) > 0
    denovolow = hl.agg.count_where(((trio_mt.proband_entry['GT'].is_het()) & trio_mt.father_entry['GT'].is_hom_ref() & trio_mt.mother_entry['GT'].is_hom_ref())) > 0 #crude de novo search w/o de novo caller. find all cases of parents hom-ref and child het. test use only. use de-novo as mode not de-novo-low for denovo calling

    xlinked = (trio_mt.locus.contig == 'chrX') & X_linked_conditions
    deNovos = filter_deNovos(deNovo, famIDarg)
    trio_mt = trio_mt.annotate_rows(dominant = hl.if_else(dominant, True, False), homrecessive = hl.if_else(homrecessive, True, False),xlinked = hl.if_else(xlinked, True, False),denovolow = hl.if_else(denovolow, True, False))
    trio_mt = trio_mt.annotate_rows(denovo = hl.if_else(hl.is_defined(deNovos[trio_mt.locus, trio_mt.alleles]), True, False))
    trio_mt = trio_mt.annotate_rows(keep_row = False)
    keep_row = trio_mt.keep_row
    for mode in inherModes:
        if mode == 'dominant':
            keep_row = hl.if_else(trio_mt.dominant == True, True, keep_row)
        if mode == 'hom-recessive':
            keep_row = hl.if_else(trio_mt.homrecessive == True, True, keep_row)
        if mode == 'x-linked':
            keep_row = hl.if_else(trio_mt.xlinked == True, True, keep_row)
        #if mode == 'samocha-de-novo-low-conf':
            #keep_row = hl.if_else(trio_mt.deNovo_low == True, True, keep_row)
        #if mode == 'samocha-de-novo-medium-conf':
            #keep_row = hl.if_else(trio_mt.deNovo_med == True, True, keep_row)
        #if mode == 'samocha-de-novo-high-conf':
            #keep_row = hl.if_else(trio_mt.deNovo_high == True, True, keep_row)
        if mode == 'de-novo':
            keep_row = hl.if_else(trio_mt.denovo == True, True, keep_row)
        if mode == 'de-novo-low':
            keep_row = hl.if_else(trio_mt.denovolow == True, True, keep_row)
    trio_mt = trio_mt.annotate_rows(keep_row = keep_row)
    trio_mt = trio_mt.filter_rows(trio_mt.keep_row == True)
    return trio_mt
##### deNovo filtering: this returns a table. Depends on custom multi allelic split. If original .vcf can be loaded /w columns in proper formats, native hail allelic split func can be used 


def filter_deNovos(deNovo, famIDarg):
    to_keep = hl.literal(famIDarg)
    deNovo = deNovo.filter(to_keep.contains(deNovo.proband.peds.famID))
    deNovo = deNovo.key_by('locus', 'alleles')
    return deNovo

def apply_gerp_score(mt, gerp_arg):
    # http://mendel.stanford.edu/SidowLab/downloads/gerp/Readme.txt Neutral rate used
    mt = mt.filter_rows(mt.gerp_scores.S > gerp_arg) 
    #mt = mt.filter_rows(mt.gerp_elements.head().S > gerp_arg).count()
    return mt
### INTERVALS
def apply_interval_filter(mt, intervals_arg):
    interval_table = None
    # get reference_genome from mt
    loci = mt.head(1).locus.collect()
    if len(loci) == 0:
        raise ValueError("input MatrixTable is emtpy")
    reference_genome = loci[0].reference_genome

    # check if intervals_arg is a file path or a list of intervals
    if os.path.isfile(intervals_arg):
        if intervals_arg.endswith('bed'):
            interval_table = hl.import_bed(intervals_arg, reference_genome=reference_genome)
        else:
            interval_table = hl.import_locus_intervals(intervals_arg, reference_genome=reference_genome)
    elif intervals_arg.startswith("gs://"):
        if intervals_arg.endswith('bed'):
            interval_table = hl.import_bed(intervals_arg, reference_genome=reference_genome)
        else:
            interval_table = hl.import_locus_intervals(intervals_arg, reference_genome=reference_genome)
    else:
        intervals_string = intervals_arg

    # parse intervals
    if interval_table:
        mt = mt.filter_rows(hl.is_defined(interval_table[mt.locus]))
       # mt = mt.annotate_rows(scATAC_source = interval_table[mt.locus].target) #annotate matrix table
    else:
        intervals = []
        mt = mt.annotate_rows(scATAC_source = '') #annotate matrix table
        for match in INTERVAL_REGEXP.finditer(intervals_string):
            interval_string = match.group(0)
            try:
                interval = hl.parse_locus_interval(interval_string, reference_genome=reference_genome)
                intervals.append(interval)
            except Exception as e:
                #logger.warning(f"Unable to parse interval {interval_string}: {e}" )
                print('error')

        if len(intervals) == 0:
            raise ValueError(f"Unable to parse any intervals from --intervals arg: {intervals_arg}")

    #logger.info(f"Parsed {len(intervals)} intervals from --intervals arg: {intervals_arg}")
    
        mt = hl.filter_intervals(mt, intervals, keep=True)
    return mt
def reformat_schema_for_output_to_tsv(mt, famFile):


    info = mt.info
    info = info.annotate(GerpScore = mt.gerp_scores['S'])
    #info = info.annotate(ClinSig = mt.clinvar_variant_summary['ClinicalSignificance']) # If using clinvar or vep annotations, uncomment these
    #info = info.annotate(TranscriptConsequences = mt.most_severe_consequence)
    info = info.annotate(topmedAF = mt.topmed_freq)
    info = info.annotate(gnomadGenomesAF = mt.gnomad_genomes_freq)
    info = info.annotate(gnomadExomesAF = mt.gnomad_exomes_freq)
    info = info.annotate(gnomadGenomesAC = mt.gnomad_genomes_AC)
    info = info.annotate(gnomadExomesAC = mt.gnomad_exomes_AC)
    info = info.select('GerpScore', 'topmedAF', 'gnomadGenomesAF', 'gnomadExomesAF', 'gnomadGenomesAC', 'gnomadExomesAC')
    mt = mt.annotate_rows(info = info)
    pedigree = hl.Pedigree.read(famFile)
    trios = hl.trio_matrix(mt, pedigree, complete_trios = False)
    trios = trios.annotate_cols(affected = trios.proband.peds.affected)
    trios = trios.annotate_cols(famID =  trios.proband.peds.famID)
    trios = trios.annotate_cols(pheno = trios.proband.pheno)
    trios = trios.key_cols_by('id', 'famID', 'affected', 'pheno')
    trios = trios.key_rows_by('locus', 'alleles', 'info') #, 'scATAC_source')
    trios_table = trios.entries().select('proband_entry', 'father_entry', 'mother_entry')
    trios_table = trios_table.filter(trios_table.proband_entry['GT'].is_non_ref()) # only return affected individuals who posses variant 
    trios_table = trios_table.annotate(GerpScore = trios_table.info.GerpScore)
    #trios_table = trios_table.annotate(ClinSig = trios_table.info.ClinSig)
    trios_table = trios_table.annotate(TopMedAF = trios_table.info.topmedAF)
    trios_table = trios_table.annotate(gnomadGenomesAF = trios_table.info.gnomadGenomesAF)
    trios_table = trios_table.annotate(gnomadExomesAF = trios_table.info.gnomadExomesAF)
    trios_table = trios_table.annotate(gnomadGenomesAC = trios_table.info.gnomadGenomesAC)
    trios_table = trios_table.annotate(gnomadExomesAC = trios_table.info.gnomadExomesAC)
    return  trios_table

def main():
    args = parse_args()
    print(args.family_id)
    mt = hl.read_matrix_table(args.input_mt_path)
    topmed = hl.read_table('topmed.ht')
    annotations = hl.read_table(args.hail_db)
    if args.intervals:
        mt = apply_interval_filter(mt, args.intervals)
        print('Intervals')

    mt.describe()
    mt = custom_split_multi_hts(mt, keep_star = False)
    mt = mt.distinct_by_row()
    #Annotate with pedigree info
    if args.fam_file:
        ped = hl.import_table(args.fam_file, no_header=True)
        ped = ped.select(famID = ped.f0, proband_ID = ped.f1, father_ID = ped.f2, mother_id= ped.f3, sex = ped.f4, affected=ped.f5)
        ped = ped.key_by('proband_ID')
        mt = mt.annotate_cols(peds=ped[mt.s]) # annotate with pedigree 
    pheno = hl.import_table('engle_ccdd_032619.ped', no_header = False)
    pheno = pheno.key_by('SAMPLE')
    mt = mt.annotate_cols(pheno=pheno[mt.s].PHENOTYPE) # annotate with pedigree 
    mt = mt.annotate_entries(AB = hl.min(mt.AD[0]/hl.sum(mt.AD), mt.AD[1]/hl.sum(mt.AD))) # AD: [Ref, Alt1, Alt2, etc]
    mt = mt.annotate_rows(topmed_freq = topmed[mt.locus, mt.alleles].info['AF'][0])
    mt = mt.annotate_rows(gnomad_genomes_freq = annotations[mt.locus, mt.alleles].gnomad_genomes_freq)
    mt = mt.annotate_rows(gnomad_exomes_freq = annotations[mt.locus, mt.alleles].gnomad_exomes_freq)
    mt = mt.annotate_rows(gnomad_genomes_AC = annotations[mt.locus, mt.alleles].gnomad_genomes_AC)
    mt = mt.annotate_rows(gnomad_exomes_AC = annotations[mt.locus, mt.alleles].gnomad_exomes_AC)
    mt = mt.annotate_rows(variant_qc = annotations[mt.locus, mt.alleles].variant_qc)
    mt = mt.annotate_rows(most_severe_consequence = annotations[mt.locus, mt.alleles].most_severe_consequence)
    mt = mt.annotate_rows(gerp_scores = annotations[mt.locus, mt.alleles].gerp_scores)
    mt = mt.annotate_rows(clinvar_variant_summary = annotations[mt.locus, mt.alleles].clinvar_variant_summary)
    mt = mt.annotate_rows(clinvar_variant_summary = annotations[mt.locus, mt.alleles].clinvar_variant_summary)
    print('annotations complete.')



    if args.pass_filter:
        mt = pass_filter(mt)
        print('pass filter')

    if args.gnomad_exomes_af:
        mt = apply_max_gnomad_exomes_AF(mt,args.gnomad_exomes_af)
        print('gnomad_exomes_AF')

    if args.gnomad_genomes_af:
        mt = apply_max_gnomad_genomes_AF(mt,args.gnomad_genomes_af)
        print('gnomad_genomes_AF')

    if args.topmed_af:
        mt = apply_topmed_AF(mt,args.topmed_af)
        print('topmed_AF')
 
    if args.gnomad_exomes_ac:
        mt = apply_max_gnomad_exomes_AC(mt, args.gnomad_exomes_ac)
        print('gnomad_exomes AC')

    if args.gnomad_genomes_ac:
        mt = apply_max_gnomad_genomes_AC(mt, args.gnomad_genomes_ac)
        print('gnomad_genomes_AC')

    if args.gq:
        mt = apply_min_GQ(mt, args.gq)
        print('min GQ')

    if args.ab:
        mt = apply_AB_filter(mt, args.ab)
        print('min AB')

    if args.gt:
        mt = apply_sample_genotype(mt, args.gt)
        print('Sample genotype')

    if args.annotations:
        mt = apply_annotations_filter(mt, args.annotations)
        print('Annotations')

    if args.gerpscore is not None:
        mt = apply_gerp_score(mt, args.gerpscore)
        print('GERP score')

    if args.genes:
        mt = filter_by_genes(mt, args.genes)
        print('Genes')

    
    if args.inheritance and args.family_id and args.fam_file:
        
        
        mt = apply_inher_filter(mt, args.fam_file, args.family_id, args.inheritance)
        print('Inheritance')
   
    

    #vars = mt.count()[0]
    #print('Unique variants before filtering: ', vars)
    to_keep = hl.literal(args.family_id)
    mt = mt.annotate_rows(unaffected_hom_var = hl.agg.count_where(((mt.peds.affected == '1') & mt.GT.is_hom_var())))
    mt = mt.annotate_rows(num_unaffected = hl.agg.count_where(((mt.peds.affected == '1') & mt.GT.is_non_ref())))

    if args.filter_unaffected:
        mt = mt.filter_rows(mt.num_unaffected == 0) #Get rid of any variants where unaffected (& not potential incomplete penetrance) have the variant; domniant searches
        #print('Unique variants after filtering: ', mt.count()[0])
    if args.filter_rec:
        mt = mt.filter_rows(mt.unaffected_hom_var > 0, keep = False) # same as filter_unaffected, but for recessive searches
        #print('Unique variants after filtering: ', mt.count()[0])

    if args.write_file: #Writing final matrix of variants
        print('storing matrix table')
        mt.write(args.write_file, overwrite = True) #immediately trim the size of the MT down
        print('reading matrix table')
        mt = hl.read_matrix_table(args.write_file)
        mt = mt.repartition(n_partitions = 50)
        print('matrix table read')
    #must be hail table to export to tsv
    print('Reformating schema to table')
    trios_table = reformat_schema_for_output_to_tsv(mt, args.fam_file) #Return individuals with variant



    if args.output_tsv:
        output_path = args.output_tsv
        #if not output_path.endswith("gz"):
            #output_path += ".gz"
    else:
       	output_path = f'{args.input_mt_path.replace(".mt", "")}.tsv.gz'

    print('Writing variants table')
    trios_table.head(args.max_variants).export(output_path, header=True)
    #mt.write('output_mt.mt', overwrite=True)
    print('Wrote variants table to ', output_path)


if __name__ == "__main__":
    main()

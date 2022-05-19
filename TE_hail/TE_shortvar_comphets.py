import hail as hl
hl.init()
union_MT = hl.read_matrix_table("TE-shortvar-union.mt")


DRS_fams = ['162','179','181','188','212','237','240','241','245','1', '6', '7', '8', '12', '17', '18', '19', '20', '26', '30', '36', '43', '51', '52', '53', '54', '55', '57', '61', '62', '65', '67', '93', '98', '100', '102', '101', '126', '130', '131', '135', '138', '143', '145', '149', '150', '151', '152', '153', '161', '166', '169', '170', '172', '189', '190', '191', '193', '211', '213', '220', '221', '222', '223', '228', '233', '235', '238', '246', '257', '259', '263', '264', '90']

moebius_fams = ['5', '31', '68', '70', '72', '73', '82', '83', '84', '85', '88', '89', '104', '105', '110', '111', '112', '113', '114', '115', '116', '117', '118', '119', '120', '121', '122', '123', '141', '142', '144', '148', '154', '157', '159', '167', '171', '174', '182', '183', '184', '185', '186', '199', '200', '201', '202', '203', '204', '206', '207', '208', '249']

CFEOM_fams = ['10', '109', '125', '128', '13', '137', '139', '147', '155', '160', '165', '197', '198', '209', '21', '22', '226', '23', '234', '244', '244', '252', '254', '255', '256', '262', '27', '33', '37', '38', '4', '41', '49', '56', '58', '71', '76', '77', '78', '79', '81', '91', '96', '99', '270', '187','239']

FNP_fams = ['2', '3', '9', '11', '28', '29', '39', '40', '42', '44', '45', '46', '50', '69']

MGJW_fams = ['268', '258', '242', '229', '195', '173', '163', '140', '136', '133', '129', '124', '108', '107', '101', '60', '47', '34', '24', '15', '14','66','95','97','146','164','176']

#no ptosis potential recessive families

CFP_fams = ['64','158','168','180','192','205','87','215']

def write_compound_hets(mt, bed, fams, output):
    ped = hl.import_table('/home/adigioia/lauren/hail/tutorials/engleccdd-2.fam', no_header=True)
    ped = ped.select(famID = ped.f0, proband_ID = ped.f1, father_ID = ped.f2, mother_id= ped.f3, sex = ped.f4, affected=ped.f5)
    ped = ped.key_by('proband_ID')
    mt = mt.annotate_cols(peds=ped[mt.s]) 
    ht = hl.import_locus_intervals(bed, reference_genome='GRCh38')
    mt = apply_interval_filter(mt, bed)
    mt = mt.annotate_rows(peak_label=ht[mt.locus].target)
    to_keep = hl.literal(fams)
    mt = mt.filter_cols(to_keep.contains(mt.peds.famID))
    pedigree = hl.Pedigree.read('/home/adigioia/lauren/hail/tutorials/engleccdd-2.fam')
    triomt = hl.trio_matrix(mt, pedigree, complete_trios=True)

    peak_mt = triomt.group_rows_by(triomt.peak_label).aggregate(
        father_proband_short = hl.agg.count_where((triomt.father_entry['GT'].is_non_ref() & triomt.proband_entry['GT'].is_het())),
        mother_proband_short =  hl.agg.count_where((triomt.mother_entry['GT'].is_non_ref() & triomt.proband_entry['GT'].is_het())),
        father_proband_TE = hl.agg.count_where((hl.is_defined(triomt.proband_entry['x']) & hl.is_defined(triomt.father_entry['x']))),
        mother_proband_TE =  hl.agg.count_where((hl.is_defined(triomt.proband_entry['x']) & hl.is_defined(triomt.mother_entry['x']))),
        loci = hl.agg.filter(((triomt.proband_entry['GT'].is_het() & (triomt.mother_entry['GT'].is_non_ref() | triomt.father_entry['GT'].is_non_ref())) | (hl.is_defined(triomt.proband_entry['x']) & (hl.is_defined(triomt.father_entry['x']) | hl.is_defined(triomt.mother_entry['x'])))),hl.agg.collect(triomt.locus)),
        alleles = hl.agg.filter(((triomt.proband_entry['GT'].is_het() & (triomt.mother_entry['GT'].is_non_ref() | triomt.father_entry['GT'].is_non_ref())) | (hl.is_defined(triomt.proband_entry['x']) & (hl.is_defined(triomt.father_entry['x']) | hl.is_defined(triomt.mother_entry['x'])))),hl.agg.collect(triomt.alleles)),
   
    )
    peak_mt = peak_mt.filter_entries(((peak_mt.father_proband_short>0) | (peak_mt.father_proband_TE>0)) & ((peak_mt.mother_proband_short> 0)|(peak_mt.mother_proband_TE>0))) #find peaks where proband inherits variant from father and mother
    peak_mt = peak_mt.filter_entries(peak_mt.alleles.length() > 1) #Eliminate case where variant from mother and father is the same variant
    peak_mt.entries().export(output)


write_compound_hets(union_MT, "CN34-intervals-padded.bed", CFEOM_fams, 'CFEOM-sv-TE-comphets.tsv')
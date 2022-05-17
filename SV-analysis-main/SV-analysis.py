#!/usr/bin/env python

## Intersections with CN34 peaks
## SV's present in unaffected individuals removed (remove SV ID)
## Results for large pedigrees filtered (remove hits in large pedigrees)
## Remove results that are present in a family w/o matching inher pattern (remove SV ID); identify variants present only 1x in family, or 2x in family where one is affected="3"
## https://github.com/lorinairs/SV-analysis/blob/main/SV-analysis.py


import argparse
import hail as hl
hl.init(default_reference='GRCh38')
def parse_args():
    p = argparse.ArgumentParser(description="aggregate filtered SV results")

    p.add_argument("-p", "--pheno", help="")
    p.add_argument("-c", "--cranialnerve", help="")

    args = p.parse_args()
    return args

args = parse_args()

infile = ''.join(['SV-results/',args.cranialnerve,'-SV-inher-exploded.tsv'])
print(infile)
allSV =''.join(['SV-results/',args.cranialnerve,'-SV-inher-aggregated.tsv'])
multiSV = ''.join(['SV-results/',args.cranialnerve,'-SV-inher-multihit.tsv'])
ht = hl.import_table(infile)
ped = hl.import_table('../engleccdd-2.fam',no_header=True)
ped = ped.rename({'f0' : 'row', 'f1' : 'samples', 'f2':'father', 'f3':'mother', 'f4':'sex', 'f5':'affected'})
ped= ped.key_by(ped.samples)

ht = ht.annotate(affected = ped[ht.samples].affected) ## annotate w affected status


### Remove SV's that don't fit custom inheritance structures
### Import SV VCF as MT
mt = hl.read_matrix_table('inhsubsetSV.mt')

valid = hl.literal(ht.aggregate(hl.agg.collect(ht.name))) # Collect all valid SV names - ht has already filtered by region w/ Genomic Ranges; use these names to filter the full MT again

mt = mt.filter_rows(valid.contains(mt.rsid)) #filter mt to only SV's present in regions

## Use inhertiance structures built for short variance, TE searches to create mt of valid SV's
#################### CUSTOM STARTS HERE
#Family224
mt224 = mt

#one_het = ((mt224.s == 'S224-7') & (mt224.GT.is_hom_ref()))
#two_het = ((mt224.s == 'S224-8') & (mt224.GT.is_hom_ref()))
#three_ref = ((mt224.s == 'S224-9') & (mt224.GT.is_hom_ref()))
seven_het = ((mt224.s == 'S224-1-2') & (mt224.GT.is_het()))
eight_het = ((mt224.s == 'S224-2') & (mt224.GT.is_het()))


#case1 = hl.agg.count_where(one_het) > 0
#case2 = hl.agg.count_where(two_het) > 0
case3 = hl.agg.count_where(seven_het) > 0
case4 = hl.agg.count_where( eight_het) > 0   
#case5 = hl.agg.count_where(three_ref) > 0 

mt224 = mt224.filter_rows(case3 & case4)

#Family260
mt260 = mt

five = ((mt260.s == 'S260-5') & (mt260.GT.is_hom_ref()))
four = ((mt260.s == 'S260-4') & (mt260.GT.is_hom_ref()))
two = ((mt260.s == 'S260-2') & (mt260.GT.is_hom_ref()))
three = ((mt260.s == 'S260-3') & (mt260.GT.is_het()))
eight = ((mt260.s == 'S260-8') & (mt260.GT.is_het()))
one = ((mt260.s == 'S260-1-2') & (mt260.GT.is_het()))


case1 = hl.agg.count_where(five) > 0
case2 = hl.agg.count_where(four) > 0
case3 = hl.agg.count_where(two) > 0
case4 = hl.agg.count_where(three) > 0   
case5 = hl.agg.count_where(eight) > 0 
case6 = hl.agg.count_where(one) > 0 

mt260 = mt260.filter_rows(case1 & case2 & case3 & case4 & case5 & case6)




#Family225
mt225 = mt

five = ((mt225.s == 'S225-5') & (mt225.GT.is_hom_ref()))
one = ((mt225.s == 'S225-1') & (mt225.GT.is_het()))
two = ((mt225.s == 'S225-2') & (mt225.GT.is_het()))
six = ((mt225.s == 'S225-6') & (mt225.GT.is_het()))


case2 = hl.agg.count_where(five) > 0
case3 = hl.agg.count_where(one) > 0
case4 = hl.agg.count_where(two) > 0   
case5 = hl.agg.count_where(six) > 0 

mt225 = mt225.filter_rows(case2 & case3 & case4 & case5)

#Family216
mt216 = mt

three = ((mt216.s == 'S216-3') & (mt216.GT.is_het()))
one = ((mt216.s == 'S216-1') & (mt216.GT.is_het()))
two = ((mt216.s == 'S216-2') & (mt216.GT.is_het()))

case1 = hl.agg.count_where(three) > 0
case2 = hl.agg.count_where(one) > 0
case3 = hl.agg.count_where(two) > 0   

mt216 = mt216.filter_rows(case1 & case2 & case3 )


#Family25
mt25 = mt 

two = ((mt25.s == 'S25-2') & (mt25.GT.is_non_ref()))
one = ((mt25.s == 'S25-1-3') & (mt25.GT.is_het()))
four = ((mt25.s == 'S25-4') & (mt25.GT.is_hom_ref()))
three = ((mt25.s == 'S25-3-2') & (mt25.GT.is_hom_ref()))



case1 = hl.agg.count_where(two) > 0
case2 = hl.agg.count_where(one) > 0
case3 = hl.agg.count_where(four) > 0
case4 = hl.agg.count_where(three) > 0   

mt25 = mt25.filter_rows(case1 & case2 & case3 & case4)


# ## CFP CUSTOM DOMINANT

#Family194
mt194 = mt 

four = ((mt194.s == 'S194-4') & (mt194.GT.is_non_ref()))
two = ((mt194.s == 'S194-2R') & (mt194.GT.is_non_ref()))
twelve = ((mt194.s == 'S194-12') & (mt194.GT.is_non_ref()))
eight = ((mt194.s == 'S194-8') & (mt194.GT.is_het()))
seven = ((mt194.s == 'S194-7') & (mt194.GT.is_het()))
one = ((mt194.s == 'S194-1') & (mt194.GT.is_non_ref()))




case1 = hl.agg.count_where(four) > 0
case2 = hl.agg.count_where(two) > 0
case3 = hl.agg.count_where(twelve) > 0
case4 = hl.agg.count_where(eight) > 0   
case5 = hl.agg.count_where(seven) > 0 
case6 = hl.agg.count_where(one) > 0 


mt194 = mt194.filter_rows(case1 & case2 & case3 & case4 & case5 & case6 )


mt156 = mt

four = ((mt156.s == 'S156-4') & (mt156.GT.is_non_ref()))
two = ((mt156.s == 'S156-2') & (mt156.GT.is_non_ref()))
three = ((mt156.s == 'S156-3') & (mt156.GT.is_non_ref()))


case1 = hl.agg.count_where(four) > 0
case3 = hl.agg.count_where(two) > 0
case4 = hl.agg.count_where(three) > 0   


mt156 = mt156.filter_rows(case1 & case3 & case4 )


#Family265
mt265 = mt 

one = ((mt265.s == 'S265-1') & (mt265.GT.is_het()))
two = ((mt265.s == 'S265-2') & (mt265.GT.is_het()))
three = ((mt265.s == 'S265-3') & (mt265.GT.is_het()))

case1 = hl.agg.count_where(one) > 0
case2 = hl.agg.count_where(two) > 0
case3 = hl.agg.count_where(three) > 0

mt265 = mt265.filter_rows(case1 & case2 & case3)


# Family 132 inheritance patterns
mt132 = mt #.filter_cols(mt.peds.famID == '132')

eight = ((mt132.s == 'S132-8') & (mt132.GT.is_non_ref()))
nine = ((mt132.s == 'S132-9-2') & (mt132.GT.is_non_ref()))
six = ((mt132.s == 'S132-6') & (mt132.GT.is_non_ref()))

case3 = hl.agg.count_where(eight) > 0
case4 = hl.agg.count_where(nine) > 0
case5 = hl.agg.count_where(six) > 0

mt132 = mt132.filter_rows(case3 & case4 & case5)


# Family 219 inheritance patterns
mt219 = mt #.filter_cols(mt.peds.famID == '219')

two = ((mt219.s == 'S219-2') & (mt219.GT.is_non_ref()))
three = ((mt219.s == 'S219-3') & (mt219.GT.is_non_ref()))
one = ((mt219.s == 'S219-1') & (mt219.GT.is_non_ref()))
five = ((mt219.s == 'S219-5') & (mt219.GT.is_non_ref()))
ten = ((mt219.s == 'S219-10') & (mt219.GT.is_non_ref()))


case1 = hl.agg.count_where(three) > 0
case2 = hl.agg.count_where(two) > 0
case3 = hl.agg.count_where(one) > 0
case4 = hl.agg.count_where(five) > 0
case5 = hl.agg.count_where(ten) > 0


#mt219 = mt219.annotate_rows(case1 = hl.if_else(case1, True, False))
mt219 = mt219.filter_rows(case1 & case2 & case3 & case4 & case5)


#Family92
mt92 = mt 

eleven = ((mt92.s == 'S92-11') & (mt92.GT.is_non_ref()))
twofour = ((mt92.s == 'S92-24') & (mt92.GT.is_het()))
three = ((mt92.s == 'S92-3') & (mt92.GT.is_het()))
two = ((mt92.s == 'S92-2') & (mt92.GT.is_hom_ref()))
one = ((mt92.s == 'S92-1') & (mt92.GT.is_het()))
six = ((mt92.s == 'S92-6') & (mt92.GT.is_het()))
nine = ((mt92.s == 'S92-9R') & (mt92.GT.is_het()))


case1 = hl.agg.count_where(eleven) > 0
case2 = hl.agg.count_where(twofour) > 0
case3 = hl.agg.count_where(three) > 0
case4 = hl.agg.count_where(two) > 0
case5 = hl.agg.count_where(one) > 0
case6 = hl.agg.count_where(six) > 0
case7 = hl.agg.count_where(nine) > 0

mt92 = mt92.filter_rows(case1 & case2 & case3 & case4 & case5 & case6 & case7 )


#Family134
mt134 = mt 

six = ((mt134.s == 'S134-6') & (mt134.GT.is_hom_ref()))
three = ((mt134.s == 'S134-3') & (mt134.GT.is_het()))
four = ((mt134.s == 'S134-4') & (mt134.GT.is_het()))
five = ((mt134.s == 'S134-5') & (mt134.GT.is_het()))
one = ((mt134.s == 'S134-1-2') & (mt134.GT.is_het()))


case1 = hl.agg.count_where(six) > 0
case2 = hl.agg.count_where(three) > 0
case3 = hl.agg.count_where(four) > 0
case4 = hl.agg.count_where(five) > 0
case5 = hl.agg.count_where(one) > 0

mt134 = mt134.filter_rows(case1 & case2 & case3 & case4 & case5)


#Family175
mt175 = mt 

one = ((mt175.s == 'S175-1') & (mt175.GT.is_het()))
eight = ((mt175.s == 'S175-8') & (mt175.GT.is_het()))
seven = ((mt175.s == 'S175-7') & (mt175.GT.is_het()))
two = ((mt175.s == 'S175-2') & (mt175.GT.is_het()))
three = ((mt175.s == 'S175-3') & (mt175.GT.is_hom_ref()))

case1 = hl.agg.count_where(one) > 0
case2 = hl.agg.count_where(eight) > 0
case3 = hl.agg.count_where(seven) > 0
case4 = hl.agg.count_where(two) > 0
case5 = hl.agg.count_where(three) > 0

mt175 = mt175.filter_rows(case1 & case2 & case3 & case4 & case5)


mt175_2 = mt 

one = ((mt175_2.s == 'S175-1') & (mt175_2.GT.is_het()))
eight = ((mt175_2.s == 'S175-8') & (mt175_2.GT.is_het()))
seven = ((mt175_2.s == 'S175-7') & (mt175_2.GT.is_hom_ref()))
two = ((mt175_2.s == 'S175-2') & (mt175_2.GT.is_hom_ref()))
three = ((mt175_2.s == 'S175-3') & (mt175_2.GT.is_hom_ref()))

case1 = hl.agg.count_where(one) > 0
case2 = hl.agg.count_where(eight) > 0
case3 = hl.agg.count_where(seven) > 0
case4 = hl.agg.count_where(two) > 0
case5 = hl.agg.count_where(three) > 0

mt175_2 = mt175_2.filter_rows(case1 & case2 & case3 & case4 & case5)


#Family178
mt178 = mt 

three = ((mt178.s == 'S178-3') & (mt178.GT.is_het()))
four = ((mt178.s == 'S178-4') & (mt178.GT.is_het()))
five = ((mt178.s == 'S178-5') & (mt178.GT.is_het()))
one = ((mt178.s == 'S178-1') & (mt178.GT.is_het()))
two = ((mt178.s == 'S178-2') & (mt178.GT.is_hom_ref()))

case1 = hl.agg.count_where(three) > 0
case2 = hl.agg.count_where(four) > 0
case3 = hl.agg.count_where(five) > 0
case4 = hl.agg.count_where(one) > 0
case5 = hl.agg.count_where(two) > 0

mt178 = mt178.filter_rows(case1 & case2 & case3 & case4 & case5)


#Family181
mt181 = mt 

one = ((mt181.s == 'S181-1') & (mt181.GT.is_het()))
three = ((mt181.s == 'S181-3') & (mt181.GT.is_het()))
two = ((mt181.s == 'S181-2') & (mt181.GT.is_hom_ref()))
four = ((mt181.s == 'S181-4') & (mt181.GT.is_hom_ref()))
five = ((mt181.s == 'S181-5') & (mt181.GT.is_hom_ref()))

case1 = hl.agg.count_where(one) > 0
case2 = hl.agg.count_where(three) > 0
case3 = hl.agg.count_where(two) > 0
case4 = hl.agg.count_where(four) > 0
case5 = hl.agg.count_where(five) > 0

mt181 = mt181.filter_rows(case1 & case2 & case3 & case4 & case5)


#Family214
mt214 = mt 

two = ((mt214.s == 'S214-2') & (mt214.GT.is_het()))
three = ((mt214.s == 'S214-3-2') & (mt214.GT.is_het()))
two = ((mt214.s == 'S214-2') & (mt214.GT.is_hom_ref()))
four = ((mt214.s == 'S214-4') & (mt214.GT.is_hom_ref()))


case1 = hl.agg.count_where(one) > 0
case2 = hl.agg.count_where(three) > 0
case3 = hl.agg.count_where(two) > 0
case4 = hl.agg.count_where(four) > 0


mt214 = mt214.filter_rows(case1 & case2 & case3 & case4 )


#Family232
mt232 = mt 

two = ((mt232.s == 'S232-2') & (mt232.GT.is_het()))
four = ((mt232.s == 'S232-4') & (mt232.GT.is_het()))
nine = ((mt232.s == 'S232-9') & (mt232.GT.is_het()))
five = ((mt232.s == 'S232-5') & (mt232.GT.is_het()))
one = ((mt232.s == 'S232-1') & (mt232.GT.is_hom_ref()))
six = ((mt232.s == 'S232-6') & (mt232.GT.is_hom_ref()))


case1 = hl.agg.count_where(two) > 0
case2 = hl.agg.count_where(four) > 0
case3 = hl.agg.count_where(nine) > 0
case4 = hl.agg.count_where(five) > 0
case5 = hl.agg.count_where(one) > 0
case6 = hl.agg.count_where(six) > 0

mt232 = mt232.filter_rows(case1 & case2 & case3 & case4 & case5 & case6)


#Family80
mt80 = mt 

one = ((mt80.s == 'S80-1-2') & (mt80.GT.is_het()))
four = ((mt80.s == 'S80-4') & (mt80.GT.is_het()))
eight = ((mt80.s == 'S80-8-2') & (mt80.GT.is_het()))
ten = ((mt80.s == 'S80-10') & (mt80.GT.is_het()))


case1 = hl.agg.count_where(one) > 0
case2 = hl.agg.count_where(four) > 0
case3 = hl.agg.count_where(eight) > 0
case4 = hl.agg.count_where(ten) > 0

mt80 = mt80.filter_rows(case1 & case2 & case3 & case4)


# ## MGJW CUSTOM DOMINANT

# Family66
mt66 = mt 

two = ((mt66.s == 'S66-2') & (mt66.GT.is_het()))
one = ((mt66.s == 'S66-1') & (mt66.GT.is_het()))
four = ((mt66.s == 'S66-4') & (mt66.GT.is_het()))


case1 = hl.agg.count_where(two) > 0
case2 = hl.agg.count_where(one) > 0
case3 = hl.agg.count_where(four) > 0

mt66 = mt66.filter_rows(case1 & case2 & case3)


# Family95
mt95 = mt 

four = ((mt95.s == 'S95-4') & (mt95.GT.is_het()))
two = ((mt95.s == 'S95-2') & (mt95.GT.is_het()))
one = ((mt95.s == 'S95-1') & (mt95.GT.is_het()))
five = ((mt95.s == 'S95-5') & (mt95.GT.is_hom_ref()))


case1 = hl.agg.count_where(four) > 0
case2 = hl.agg.count_where(two) > 0
case3 = hl.agg.count_where(one) > 0
case4 = hl.agg.count_where(five) > 0

mt95 = mt95.filter_rows(case1 & case2 & case3 & case4)


# Family97
mt97 = mt 

two = ((mt97.s == 'S97-2') & (mt97.GT.is_het()))
one = ((mt97.s == 'S97-1') & (mt97.GT.is_het()))
four = ((mt97.s == 'S97-4') & (mt97.GT.is_het()))

case1 = hl.agg.count_where(four) > 0
case2 = hl.agg.count_where(two) > 0
case3 = hl.agg.count_where(one) > 0

mt97 = mt97.filter_rows(case1 & case2 & case3)


#Family 146
mt146 = mt 

one = ((mt146.s == 'S146-1') & (mt146.GT.is_het()))
two = ((mt146.s == 'S146-2') & (mt146.GT.is_het()))
four = ((mt146.s == 'S146-4-3') & (mt146.GT.is_het()))
three = ((mt146.s == 'S146-3') & (mt146.GT.is_hom_ref()))

case1 = hl.agg.count_where(one) > 0
case2 = hl.agg.count_where(two) > 0
case3 = hl.agg.count_where(four) > 0
case4 = hl.agg.count_where(three) > 0


mt146 = mt146.filter_rows(case1 & case2 & case3 & case4)


# Family 164
mt164 = mt 

one = ((mt164.s == 'S164-1') & (mt164.GT.is_het()))
five = ((mt164.s == 'S164-5') & (mt164.GT.is_het()))
three = ((mt164.s == 'S164-3') & (mt164.GT.is_hom_ref()))
two = ((mt164.s == 'S164-2') & (mt164.GT.is_hom_ref()))
four = ((mt164.s == 'S164-4') & (mt164.GT.is_hom_ref()))

case1 = hl.agg.count_where(one) > 0
case2 = hl.agg.count_where(five) > 0
case3 = hl.agg.count_where(three) > 0
case4 = hl.agg.count_where(two) > 0
case5 = hl.agg.count_where(four) > 0

mt164 = mt164.filter_rows(case1 & case2 & case3 & case4 & case5)


# Family 176
mt176 = mt 

two = ((mt176.s == 'S176-2') & (mt176.GT.is_het()))
one = ((mt176.s == 'S176-1') & (mt176.GT.is_het()))
three = ((mt176.s == 'S176-3') & (mt176.GT.is_hom_ref()))

case1 = hl.agg.count_where(two) > 0
case2 = hl.agg.count_where(one) > 0
case3 = hl.agg.count_where(three) > 0

mt176 = mt176.filter_rows(case1 & case2 & case3 )


## Ptosis

#Family 248
mt248 = mt

one = ((mt248.s == 'S248-1-4') & (mt248.GT.is_het()))
four = ((mt248.s == 'S248-4') & (mt248.GT.is_het()))
seven = ((mt248.s == 'S248-7-2') & (mt248.GT.is_het()))
fourteen = ((mt248.s == 'S248-14-2') & (mt248.GT.is_het()))
twenty = ((mt248.s == 'S248-20-3') & (mt248.GT.is_het()))


case1 = hl.agg.count_where(one) > 0
case2 = hl.agg.count_where(four) > 0
case3 = hl.agg.count_where(seven) > 0
case4 = hl.agg.count_where(fourteen) > 0
case5 = hl.agg.count_where(twenty) > 0

mt248 = mt248.filter_rows(case1 & case2 & case3 & case4 &case5)


#Family 243
mt243 = mt #.filter_cols(mt.peds.famID == '92')

two = ((mt243.s == 'S243-2') & (mt243.GT.is_het()))
one = ((mt243.s == 'S243-1') & (mt243.GT.is_hom_ref()))
three = ((mt243.s == 'S243-3') & (mt243.GT.is_hom_ref()))
four = ((mt243.s == 'S243-4') & (mt243.GT.is_het()))


case1 = hl.agg.count_where(two) > 0
case2 = hl.agg.count_where(three) > 0
case3 = hl.agg.count_where(one) > 0
case4 = hl.agg.count_where(four) > 0


mt243 = mt243.filter_rows(case1 & case2 & case3 & case4 )


# Family 253
mt253 = mt 

two = ((mt253.s == 'S253-2') & (mt253.GT.is_het()))
five = ((mt253.s == 'S253-5') & (mt253.GT.is_het()))
four = ((mt253.s == 'S253-4') & (mt253.GT.is_het()))
three = ((mt253.s == 'S253-3') & (mt253.GT.is_het()))
one = ((mt253.s == 'S253-1') & (mt253.GT.is_hom_ref()))

case1 = hl.agg.count_where(two) > 0
case2 = hl.agg.count_where(five) > 0
case3 = hl.agg.count_where(four) > 0
case4 = hl.agg.count_where(three) > 0
case5 = hl.agg.count_where(one) > 0

mt253 = mt253.filter_rows(case1 & case2 & case3 & case4 & case5)


# Family 247
mt247 = mt 

one = ((mt247.s == 'S247-1-2') & (mt247.GT.is_het()))
five = ((mt247.s == 'S247-5-2') & (mt247.GT.is_het()))
four = ((mt247.s == 'S247-4-2') & (mt247.GT.is_het()))
six = ((mt247.s == 'S247-6-2') & (mt247.GT.is_het()))
three = ((mt247.s == 'S247-3-2') & (mt247.GT.is_hom_ref()))
two = ((mt247.s == 'S247-2') & (mt247.GT.is_hom_ref()))


case1 = hl.agg.count_where(one) > 0
case2 = hl.agg.count_where(five) > 0
case3 = hl.agg.count_where(four) > 0
case4 = hl.agg.count_where(six) > 0
case5 = hl.agg.count_where(three) > 0
case6 = hl.agg.count_where(two) > 0

mt247 = mt247.filter_rows(case1 & case2 & case3 & case4 & case5 & case6)


# Family 210
mt210 = mt 

two = ((mt210.s == 'S210-2') & (mt210.GT.is_het()))
four = ((mt210.s == 'S210-4') & (mt210.GT.is_het()))
three = ((mt210.s == 'S210-3') & (mt210.GT.is_hom_ref()))

case1 = hl.agg.count_where(two) > 0
case2 = hl.agg.count_where(four) > 0
case3 = hl.agg.count_where(three) > 0

mt210 = mt210.filter_rows(case1 & case2 & case3)


#Family 218
mt218 = mt 

seven = ((mt218.s == 'S218-7') & (mt218.GT.is_het()))
eight = ((mt218.s == 'S218-8') & (mt218.GT.is_het()))
three = ((mt218.s == 'S218-3_1') & (mt218.GT.is_het()))

case1 = hl.agg.count_where(seven) > 0
case2 = hl.agg.count_where(eight) > 0
case3 = hl.agg.count_where(three) > 0

mt218 = mt218.filter_rows(case1 & case2 & case3)


# Merge MT's 
merge_all = [mt25, mt216, mt224, mt225, mt260, mt194, mt265, mt156, mt92, mt132, mt134, mt175, mt178, mt181, mt214, mt219, mt232,
 mt80,mt66, mt95, mt97, mt146, mt164, mt176,mt210, mt218, mt243, mt247, mt248, mt253]

merged = hl.MatrixTable.union_rows(*merge_all)
################# CUSTOM ENDS HERE

## List of families custom pedigrees apply to
custom_fams = hl.literal(['25','216', '224', '225', '260', '194', '265', '156', '92', '132', '134','175', '178','181','214',
              '219', '232', '80','66','95','97','146','164','176','210','218','243','247','248','253'])

custom_fam_sv = hl.literal(merged.aggregate_rows(hl.agg.collect_as_set(merged.rsid))) ## List of SV names valid for custom pedigrees - found with custom structures above

all_sv_ids = hl.literal(ht.aggregate(hl.agg.collect_as_set(ht.name))) # full list of SV names

not_custom_fam_sv = all_sv_ids.difference(custom_fam_sv) # invalid SV names for any custom family 

ht= ht.filter(custom_fams.contains(ht.famID) & (not_custom_fam_sv.contains(ht.name)), keep = False) ## filter the exploded list to remove rows from custom fams containing invalid SV's

## Aggregate results by SV
SVagg = ht.group_by(ht.SVloc).aggregate(peaks = hl.agg.collect_as_set(ht.peak),
                                        samples = hl.agg.collect_as_set(ht.samples),
                                        fams = hl.agg.collect_as_set(ht.famID),
                                        pheno = hl.agg.collect_as_set(ht.pheno),
                                        affected = hl.agg.collect_as_set(ht.affected),
                                        width = hl.agg.collect(ht.width)[0],
                                        name = hl.agg.collect(ht.name)[0],
                                        svtype = hl.agg.collect(ht.svtype)[0],
                                        algo = hl.agg.collect(ht.ALGORITHMS)[0],
                                        evidence = hl.agg.collect(ht.EVIDENCE)[0],
                                        SVlen = hl.agg.collect(ht.SVLEN)[0],
                                        SVtype = hl.agg.collect(ht.SVTYPE)[0],
                                        unresolved_type = hl.agg.collect(ht.UNRESOLVED_TYPE)[0],
                                        protein_coding_lof = hl.agg.collect(ht.PROTEIN_CODING__LOF)[0],
                                        lincrna_lof = hl.agg.collect(ht.LINCRNA__LOF)[0],
                                        protein_coding_dup_lof = hl.agg.collect(ht.PROTEIN_CODING__DUP_LOF)[0],
                                        lincrna_dup_lof = hl.agg.collect(ht.LINCRNA__DUP_LOF)[0],
                                        protein_coding_copy_gain = hl.agg.collect(ht.PROTEIN_CODING__COPY_GAIN)[0],
                                        lincrna_copy_gain = hl.agg.collect(ht.LINCRNA__COPY_GAIN)[0],
                                        protein_coding_dup_partial = hl.agg.collect(ht.PROTEIN_CODING__DUP_PARTIAL)[0],
                                        lincrna_dup_partial = hl.agg.collect(ht.LINCRNA__DUP_PARTIAL)[0],
                                        protein_coding_msv_exon_ovr = hl.agg.collect(ht.PROTEIN_CODING__MSV_EXON_OVR)[0],
                                        lincrna_msv_exon_ovr = hl.agg.collect(ht.LINCRNA__MSV_EXON_OVR)[0],
                                        protein_coding_intronic = hl.agg.collect(ht.PROTEIN_CODING__INTRONIC)[0],
                                        lincrna_intronic = hl.agg.collect(ht.LINCRNA__INTRONIC)[0],
                                        protein_coding_inv_span = hl.agg.collect(ht.PROTEIN_CODING__INV_SPAN)[0],
                                        lincrna_inv_span = hl.agg.collect(ht.LINCRNA__INV_SPAN)[0],
                                        protein_coding_utr = hl.agg.collect(ht.PROTEIN_CODING__UTR)[0],
                                        lincrna_unr = hl.agg.collect(ht.LINCRNA__UTR)[0],
                                        noncoding_span = hl.agg.collect(ht.NONCODING_SPAN)[0],
                                        noncoding_breakpoint = hl.agg.collect(ht.NONCODING_BREAKPOINT)[0],
                                        protein_coding_nearest_tss = hl.agg.collect(ht.PROTEIN_CODING__NEAREST_TSS)[0],
                                        protein_coding_intergenic = hl.agg.collect(ht.PROTEIN_CODING__INTERGENIC)[0],
                                        protein_coding_promoter = hl.agg.collect(ht.PROTEIN_CODING__PROMOTER)[0],
                                        AC = hl.agg.collect(ht.AC)[0],
                                        AN = hl.agg.collect(ht.AN)[0]
                                        )


SVagg = SVagg.filter(SVagg.affected.contains("1"), keep = False) # further remove SV's present in any unaffected individuals (filter any invalid inher variants in simple dom trios)
SVagg = SVagg.filter(SVagg.samples.length() > 1, keep = True) # remove any SV's where custom pedigree variant filtering removed inher relationship, leaving singletons
filtered_SV_ids = hl.literal(SVagg.aggregate(hl.agg.collect_as_set(SVagg.name))) # full list of SV names

ht= ht.filter(filtered_SV_ids.contains(ht.name), keep = True) 

## reaggregate, by SVLoc & family
SVagg = ht.group_by(ht.SVloc, ht.famID).aggregate(peaks = hl.agg.collect_as_set(ht.peak),
                                        samples = hl.agg.collect_as_set(ht.samples),
                                        fams = hl.agg.collect_as_set(ht.famID),
                                        pheno = hl.agg.collect_as_set(ht.pheno),
                                        affected = hl.agg.collect_as_set(ht.affected),
                                        width = hl.agg.collect(ht.width)[0],
                                        name = hl.agg.collect(ht.name)[0],
                                        svtype = hl.agg.collect(ht.svtype)[0],
                                        algo = hl.agg.collect(ht.ALGORITHMS)[0],
                                        evidence = hl.agg.collect(ht.EVIDENCE)[0],
                                        SVlen = hl.agg.collect(ht.SVLEN)[0],
                                        SVtype = hl.agg.collect(ht.SVTYPE)[0],
                                        unresolved_type = hl.agg.collect(ht.UNRESOLVED_TYPE)[0],
                                        protein_coding_lof = hl.agg.collect(ht.PROTEIN_CODING__LOF)[0],
                                        lincrna_lof = hl.agg.collect(ht.LINCRNA__LOF)[0],
                                        protein_coding_dup_lof = hl.agg.collect(ht.PROTEIN_CODING__DUP_LOF)[0],
                                        lincrna_dup_lof = hl.agg.collect(ht.LINCRNA__DUP_LOF)[0],
                                        protein_coding_copy_gain = hl.agg.collect(ht.PROTEIN_CODING__COPY_GAIN)[0],
                                        lincrna_copy_gain = hl.agg.collect(ht.LINCRNA__COPY_GAIN)[0],
                                        protein_coding_dup_partial = hl.agg.collect(ht.PROTEIN_CODING__DUP_PARTIAL)[0],
                                        lincrna_dup_partial = hl.agg.collect(ht.LINCRNA__DUP_PARTIAL)[0],
                                        protein_coding_msv_exon_ovr = hl.agg.collect(ht.PROTEIN_CODING__MSV_EXON_OVR)[0],
                                        lincrna_msv_exon_ovr = hl.agg.collect(ht.LINCRNA__MSV_EXON_OVR)[0],
                                        protein_coding_intronic = hl.agg.collect(ht.PROTEIN_CODING__INTRONIC)[0],
                                        lincrna_intronic = hl.agg.collect(ht.LINCRNA__INTRONIC)[0],
                                        protein_coding_inv_span = hl.agg.collect(ht.PROTEIN_CODING__INV_SPAN)[0],
                                        lincrna_inv_span = hl.agg.collect(ht.LINCRNA__INV_SPAN)[0],
                                        protein_coding_utr = hl.agg.collect(ht.PROTEIN_CODING__UTR)[0],
                                        lincrna_unr = hl.agg.collect(ht.LINCRNA__UTR)[0],
                                        noncoding_span = hl.agg.collect(ht.NONCODING_SPAN)[0],
                                        noncoding_breakpoint = hl.agg.collect(ht.NONCODING_BREAKPOINT)[0],
                                        protein_coding_nearest_tss = hl.agg.collect(ht.PROTEIN_CODING__NEAREST_TSS)[0],
                                        protein_coding_intergenic = hl.agg.collect(ht.PROTEIN_CODING__INTERGENIC)[0],
                                        protein_coding_promoter = hl.agg.collect(ht.PROTEIN_CODING__PROMOTER)[0],
                                        AC = hl.agg.collect(ht.AC)[0],
                                        AN = hl.agg.collect(ht.AN)[0]
                                        )

filter_cond = (((SVagg.samples.length() == 2) & SVagg.affected.contains("3")) |(SVagg.samples.length() == 1))
SVagg = SVagg.annotate(trueMatch = 
hl.if_else(filter_cond, "no",  "yes"))

SVagg = SVagg.group_by(SVagg.SVloc).aggregate(trueMatch = hl.agg.collect_as_set(SVagg.trueMatch)
                                        , name = hl.agg.collect(SVagg.name))

SVagg = SVagg.filter(SVagg.trueMatch.contains('no'), keep = False)

filtered_SV_ids = hl.literal(SVagg.aggregate(hl.agg.collect_as_set(SVagg.name[0]))) # full list of SV names
ht= ht.filter(filtered_SV_ids.contains(ht.name), keep = True) 

## Final reaggregation of filtered exploded list

SVagg = ht.group_by(ht.SVloc).aggregate(peaks = hl.agg.collect_as_set(ht.peak),
                                        samples = hl.agg.collect_as_set(ht.samples),
                                        fams = hl.agg.collect_as_set(ht.famID),
                                        pheno = hl.agg.collect_as_set(ht.pheno),
                                        affected = hl.agg.collect_as_set(ht.affected),
                                        width = hl.agg.collect(ht.width)[0],
                                        name = hl.agg.collect(ht.name)[0],
                                        svtype = hl.agg.collect(ht.svtype)[0],
                                        algo = hl.agg.collect(ht.ALGORITHMS)[0],
                                        evidence = hl.agg.collect(ht.EVIDENCE)[0],
                                        SVlen = hl.agg.collect(ht.SVLEN)[0],
                                        SVtype = hl.agg.collect(ht.SVTYPE)[0],
                                        unresolved_type = hl.agg.collect(ht.UNRESOLVED_TYPE)[0],
                                        protein_coding_lof = hl.agg.collect(ht.PROTEIN_CODING__LOF)[0],
                                        lincrna_lof = hl.agg.collect(ht.LINCRNA__LOF)[0],
                                        protein_coding_dup_lof = hl.agg.collect(ht.PROTEIN_CODING__DUP_LOF)[0],
                                        lincrna_dup_lof = hl.agg.collect(ht.LINCRNA__DUP_LOF)[0],
                                        protein_coding_copy_gain = hl.agg.collect(ht.PROTEIN_CODING__COPY_GAIN)[0],
                                        lincrna_copy_gain = hl.agg.collect(ht.LINCRNA__COPY_GAIN)[0],
                                        protein_coding_dup_partial = hl.agg.collect(ht.PROTEIN_CODING__DUP_PARTIAL)[0],
                                        lincrna_dup_partial = hl.agg.collect(ht.LINCRNA__DUP_PARTIAL)[0],
                                        protein_coding_msv_exon_ovr = hl.agg.collect(ht.PROTEIN_CODING__MSV_EXON_OVR)[0],
                                        lincrna_msv_exon_ovr = hl.agg.collect(ht.LINCRNA__MSV_EXON_OVR)[0],
                                        protein_coding_intronic = hl.agg.collect(ht.PROTEIN_CODING__INTRONIC)[0],
                                        lincrna_intronic = hl.agg.collect(ht.LINCRNA__INTRONIC)[0],
                                        protein_coding_inv_span = hl.agg.collect(ht.PROTEIN_CODING__INV_SPAN)[0],
                                        lincrna_inv_span = hl.agg.collect(ht.LINCRNA__INV_SPAN)[0],
                                        protein_coding_utr = hl.agg.collect(ht.PROTEIN_CODING__UTR)[0],
                                        lincrna_unr = hl.agg.collect(ht.LINCRNA__UTR)[0],
                                        noncoding_span = hl.agg.collect(ht.NONCODING_SPAN)[0],
                                        noncoding_breakpoint = hl.agg.collect(ht.NONCODING_BREAKPOINT)[0],
                                        protein_coding_nearest_tss = hl.agg.collect(ht.PROTEIN_CODING__NEAREST_TSS)[0],
                                        protein_coding_intergenic = hl.agg.collect(ht.PROTEIN_CODING__INTERGENIC)[0],
                                        protein_coding_promoter = hl.agg.collect(ht.PROTEIN_CODING__PROMOTER)[0],
                                        AC = hl.agg.collect(ht.AC)[0],
                                        AN = hl.agg.collect(ht.AN)[0]
                                        )

# aggregate by peak
peak_agg = ht.group_by(ht.peak).aggregate(SVloc = hl.agg.collect_as_set(ht.SVloc),
                                        samples = hl.agg.collect_as_set(ht.samples),
                                        fams = hl.agg.collect_as_set(ht.famID),
                                        pheno = hl.agg.collect_as_set(ht.pheno),
                                        width = hl.agg.collect(ht.width),
                                        name = hl.agg.collect(ht.name),
                                        svtype = hl.agg.collect(ht.svtype),
                                        AC = hl.agg.collect(ht.AC)
                                        )
peak_agg = peak_agg.filter(peak_agg.fams.length() >1, keep = True) #filter peak aggregations to collect only multiple families, single phenos
peak_agg = peak_agg.filter(peak_agg.pheno.length() >1, keep = False)

SVagg.export(allSV) #write both aggregations to table
peak_agg.export(multiSV)


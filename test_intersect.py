from multiprocessing import Pool
from argparse import ArgumentParser
import tabix

def read_file(file):
    tb=tabix.open(file)
    vcf=[]
    for x,y in chromosome_lengths.items():
        record=tb.query(str(x),1,y)  
        vcf.append(list(record))
    return vcf
def findVariant(variants):
    holdme=[]
    hold_dict={}
    for f in variants:
        for thing in f:
            holdme.append(thing[2])
            hold_dict[thing[2]]=thing
    return set(holdme), hold_dict
def findVariant_parallel(variants):
    pool=Pool()
    results=pool.map(findVariant,[file for file in variants])
    return results
def process_files_parallel(files):
    pool=Pool()
    results=pool.map(read_file,[file for file in files])
    return results

def find_intersection(found_variants,name):
    hc_var=found_variants[0][0]
    dv_var=found_variants[1][0]
    dv_dict=found_variants[1][1]
    combined=list(dv_var.intersection(hc_var))

    with open(name+'_combined_out.vcf','w') as out:
        for f in combined:
            out.write("\t".join(dv_dict[f])+'\n')

parser = ArgumentParser(description='Grab info for filter')
parser.add_argument('-c','--chromosome_lengths',help='Optional tab delimited chromosome file, if not using GRCh38')
parser.add_argument('-g','--gatk',help='Path to gatk output')
parser.add_argument('-d','--deepvariant',help='Path to DeepVariant output')
parser.add_argument('-r','--ref',help='First line of reference file to see if chr is in fasta')

args = parser.parse_args()
chromosome_lengths={}
chrom=0
if 'chr' in args.ref:
    chrom=1
if args.chromosome_lengths =='None' and chrom==0:
    chromosome_lengths={1:248956422,2:242193529,3:198295559,4:190214555,5:181538259,6:170805979,7:159345973,8:145138636,9:138394717,10:133797422,11:135086622,12:133275309,13:114364328,14:107043718,15:101991189,16:90338345,17:83257441,18:80373285,19:58617616,20:64444167,21:46709983,22:50818468,'X':156040895}
elif args.chromosome_lengths =='None' and chrom==1:
    chromosome_lengths={'chr1':248956422,'chr2':242193529,'chr3':198295559,'chr4':190214555,'chr5':181538259,'chr6':170805979,'chr7':159345973,'chr8':145138636,'chr9':138394717,'chr10':133797422,'chr11':135086622,'chr12':133275309,'chr13':114364328,'chr14':107043718,'chr15':101991189,'chr16':90338345,'chr17':83257441,'chr18':80373285,'chr19':58617616,'chr20':64444167,'chr21':46709983,'chr22':50818468,'chrX':156040895}

else:

    with open(args.chromosome_lengths) as input:
        for line in input:
            data=line.strip().split('\t')
            chromosome_lengths[data[0]]=int(data[1])
name=args.gatk.split('.glnexus_denovo_actual.hc.vcf.gz')[0]
input=[args.gatk,args.deepvariant]
read_files=process_files_parallel(input)

found_variants=findVariant_parallel(read_files)
find_intersection(found_variants,name)


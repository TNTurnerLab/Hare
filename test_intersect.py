from multiprocessing import Pool
from argparse import ArgumentParser
import tabix

def read_file(file):
    tb=tabix.open(file)
    vcf=[]
    #print(chromosome_lengths)
    for x,y in chromosome_lengths.items():
        record=tb.query("chr"+str(x),1,y)  
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
    #print(dv_dict)
    combined=list(hc_var.intersection(dv_var))
    #print(len(combined))
    combined=list(dv_var.intersection(hc_var))
    #print(len(combined))
    with open(name+'_combined_out.vcf','w') as out:
        for f in combined:
            out.write("\t".join(dv_dict[f])+'\n')


parser = ArgumentParser(description='Grab info for filter')
parser.add_argument('-f','--chromosome_lengths',help='Optional tab delimited chromosome file, if not using GRCh38',default=False)
parser.add_argument('-g','--gatk',help='Path to gatk output')
parser.add_argument('-d','--deepvariant',help='Path to DeepVariant output')
args = parser.parse_args()
chromosome_lengths={}
if args.chromosome_lengths ==False:
    chromosome_lengths={1:248956422,2:242193529,3:198295559,4:190214555,5:181538259,6:170805979,7:159345973,8:145138636,9:138394717,10:133797422,11:135086622,12:133275309,13:114364328,14:107043718,15:101991189,16:90338345,17:83257441,18:80373285,19:58617616,20:64444167,21:46709983,22:50818468,'X':156040895}
else:
    with open(args.chromosome_lengths) as input:
        for line in input:
            data=line.strip().split('\t')
            if data[0].startswith('chr'):

                chromosome_lengths[int(data[0].split('chr')[1])]=int(data[1])
            else:
                chromosome_lengths[int(data[0])]=int(data[1])
name=args.gatk.split('.glnexus_denovo_actual.hc.vcf.gz')[0]
print(name)    
input=[args.gatk,args.deepvariant]
read_files=process_files_parallel(input)

found_variants=findVariant_parallel(read_files)
find_intersection(found_variants,name)

#print(len(hello))

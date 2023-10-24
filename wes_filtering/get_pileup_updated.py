from argparse import ArgumentParser
import os
from pybedtools import BedTool
def grab_coord(file,coord,ref,outdir):
    ref_hold=''
    holdme={}
    with open(ref) as input:
        for line in input:
            if line.startswith('>chr'):
                ref_hold='chr'
                break
    print(ref_hold)
    with open(outdir+'/'+file+'_bed_coor.txt','w') as out:
        for thing in coord:
           
            if not thing.startswith('#'):
                data=thing.strip().split('\t')
                pos=data[1]
                start=str(int(data[1])-1)
                out.write(ref_hold+data[0]+'\t'+start+'\t'+pos+'\n')
                if len(data[3])==1 and len(data[4])==1:
                    holdme[ref_hold+data[0]+'_'+pos]=[data[4],thing.strip()]
                else:
                    if len(data[3]) > len(data[4]):
                        holdme[ref_hold+data[0]+'_'+pos]=[data[3][1:],thing.strip()]
                  
                    elif len(data[3]) < len(data[4]):
                        holdme[ref_hold+data[0]+'_'+pos]=[data[4][1:],thing.strip()]
                
    return(holdme)

def grab_family(famloc,family):
    pathways=''
    for f in family:
        pathways+=famloc[f]+'\t'
    return(pathways)
def read_locations(file):
    holdme={}
    with open(file) as input:
        for line in input:
            data=line.strip().split('/')
            holdme[data[-1].split('.')[0]]=line.strip()
    return(holdme)

def specific_check(father,mother,checkme,lower,upper,case_specific):
    if case_specific==False:
        check_total=lower+upper        
        if father.lower().count(checkme.lower())<=check_total and mother.lower().count(checkme.lower())<=check_total:
            return(True)
        else:
            return(False)
    if case_specific==True:

        if father.count(checkme.lower())<=lower and mother.count(checkme.lower())<=lower and father.count(checkme.upper())<=upper and mother.count(checkme.upper())<=upper:
            return(True)
        else:
            return(False)

def run_mpileup(bed,ref,family,pathways,vcffile,outdir,case_specific,lower,upper,interval):
    coord=[]
    header=''
    with open(vcffile) as input:
        for line in input:
      
            if line.startswith('#'):
                header+=line
            else:
                coord.append(line)
    kid=family[-1]
    cram_files=grab_family(pathways,family)

    check_alt=grab_coord(family[-1],coord,ref,outdir)

    os.system('/opt/conda/bin/samtools mpileup --fasta-ref '+ref+' --positions '+bed +' --output '+outdir+'/'+kid+'_mpileup '+cram_files)

    with open(outdir+'/'+kid+'_mpileup') as input:
        with open(outdir+'/'+kid+'_mpileup_filtered_checkme.txt','w') as out:
       
            hold_high,hold_low=findIntersect(outdir+'/'+family[-1]+'_bed_coor.txt',interval)
       
            for line in input:
                data=line.strip().split('\t')
                chr=data[0].split('chr')[1]
      
                checkme=check_alt[data[0]+'_'+data[1]][0]
     
                confidence=getCon(data[0]+'_'+data[1],hold_high,hold_low)
          
                if checkme.lower() in data[10].lower():
                    if case_specific==False and lower==0 and upper==0:
                        if checkme.lower() not in data[4].lower() and checkme.lower() not in data[7].lower():
                            
                            out.write(kid+'\t'+check_alt[data[0]+'_'+data[1]][1].split('\t')[2]+'\t'+'Pass\t'+confidence+'\n')
                        else:
                            out.write(kid+'\t'+check_alt[data[0]+'_'+data[1]][1].split('\t')[2]+'\t'+'Fail\t'+confidence+'\n')
                    else:
                        if specific_check(data[4],data[7],checkme,lower,upper,case_specific)==True:
                            out.write(kid+'\t'+check_alt[data[0]+'_'+data[1]][1].split('\t')[2]+'\t'+'Pass\t'+confidence+'\n')
                        else:
                            out.write(kid+'\t'+check_alt[data[0]+'_'+data[1]][1].split('\t')[2]+'\t'+'Fail\t'+confidence+'\n')
                else:
                    out.write(kid+'\t'+check_alt[data[0]+'_'+data[1]][1].split('\t')[2]+'\t'+'Fail\t'+confidence+'\n')
def getCon(checkme,high,low):
    if checkme in high:
        return('high_confidence')
    elif checkme in low:
        return('low_confidence')
def checkCon(compare):
    holdme=set()               
    for i in range(len(compare)):
        end=str(compare[i].stop)
        chrom=compare[i].chrom
        holdme.add(chrom+'_'+end)
    return(holdme)
def findIntersect(file,compare):
    a = BedTool(file)

    hold_highcon=checkCon(a.intersect(compare))
    hold_lowcon=checkCon(a.subtract(compare))
    return(hold_highcon,hold_lowcon)

  

parser = ArgumentParser(description='Grab info for filter')
parser.add_argument('-v','--vcffile',help='VCF HAT output file')
parser.add_argument('-p','--cramlocations',help='cram location')
parser.add_argument('-f','--family',help='Comma delimited family in father, mother, kid order')
parser.add_argument('-r','--reference',help='Reference file for cram/bam')
parser.add_argument('-o','--outdir',help='path to output')
parser.add_argument('-i','--interval',help='path to interval file')
parser.add_argument('-l','--lower',help='Details level of scurnity for filter. Please denote how many occurances of the alt allele at the lower case value is allowed when filtering.  Default is 0 and is not case specific',default=0,type=int)
parser.add_argument('-u','--upper',help='Details level of scurnity for filter. Please denote how many occurances of the alt allele at the upper case value is allowed when filtering.  Default is 0 and is not case specific',default=0,type=int)
parser.add_argument('-c','--case_specific',help='When filtering for alt, check to see if it is upper (high quality), or lower (low quailty) case read.  Default is not case specific.  Please ensure that you use this flag when specifying multiple quality types',action='store_true')
args = parser.parse_args()
if args.lower != 0 and args.upper !=0 and args.case_specific==False:
    print('WARNING:  If you want to define individual filter limits for high and low quality reads, --case_specific flag must be used')
    print('This run will continue under the assumption that it is NOT quality specific')

family=args.family.split(',')
fam_locations=read_locations(args.cramlocations)
print(args.case_specific)
run_mpileup(args.outdir+'/'+family[-1]+'_bed_coor.txt',args.reference,family,fam_locations,args.vcffile,args.outdir,args.case_specific,args.lower,args.upper,args.interval)
    


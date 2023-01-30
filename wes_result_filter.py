from argparse import ArgumentParser
from pybedtools import BedTool
import os

#Modifies the interval file by shrinking the positions by user input, default to 40bp.
def modifyIntervalFile(file,change):
    with open(file+'_pm'+str(change)+'_changed.bed','w') as out:
        with open(file) as input:
            for line in input:
                data=line.strip().split('\t')
                if int(data[2])-change > int(data[1])+change:
                    out.write(data[0]+'\t'+str(int(data[1])+change)+'\t'+str(int(data[2])-change)+'\n')

#Grabs output from HAT, either from direct output or aggregated folder of vcf files. 
def grabOutput(path):
    with open('combined_dnv.bed','w') as out:
        files=os.listdir(path)
        
        for file in files:
            checkme=0
            if os.path.isdir(path+'/'+file):
                vcfpath=path+'/'+file+'/'+file+'.glnexus.family.combined_intersection_filtered_gq_20_depth_10_position.vcf' 
            elif file.endswith('_position.vcf'):
                vcfpath=path+'/'+file
            else:
                vcfpath=path+'/'+file
                checkme=1
            if os.path.exists(vcfpath) and checkme==0:

                with open(vcfpath) as input:
                    
                    sample=vcfpath.split('/')[-1].split('.')[0]
                    for line in input:
                        if not line.startswith("#"):
                            data=line.strip().split('\t')
                            out.write(data[0]+'\t'+str(int(data[1])-1)+'\t'+data[1]+'\t'+sample+'_'+data[2]+'\n')
            else:
                print('The file '+vcfpath+' does not exist or is not the correct HAT output file. This file will be skipped')
def writeFullOutput(file):
    name=file.split('.bed')[0]
    with open(name+'_table.txt','w') as out:
        out.write('SAMPLE\tCHROM\tPOS_B38\tREFERENCE\tALTERNATE\tID\n')
        with open(file) as input:
            for line in input:
                data=line.strip().split('\t')
                sample=data[-1].split('_')[0]
                ref=data[-1].split('_')[3]
                alt=data[-1].split('_')[4]
                id='_'.join(data[-1].split('_')[1:])
                out.write(sample+'\t'+data[0]+'\t'+data[2]+'\t'+ref+'\t'+alt+'\t'+id+'\n')
#Calls pyBedtools to do the intersection and subtraction calls
def findIntersect(file,change,outpath):
    a = BedTool('combined_dnv.bed')
    a.intersect(file+'_pm'+str(change)+'_changed.bed').saveas(outpath+'high_confidence_dnvs.bed')
    a.subtract(file+'_pm'+str(change)+'_changed.bed').saveas(outpath+'low_confidence_dnvs.bed')
    writeFullOutput(outpath+'high_confidence_dnvs.bed')
    writeFullOutput(outpath+'low_confidence_dnvs.bed')
parser = ArgumentParser(description='Grab info for filter')
parser.add_argument('-p','--path',help='Path to out folder. This can be either a the output folder from HAT or a singular directory containing the vcf files from HAT with suffix _position.vcf',required=True)
parser.add_argument('-i','--interval',help='Path to interval file',required=True)
parser.add_argument('-n','--number',help='Number of bp to change the interval file positions by',default=40,type=int)
parser.add_argument('-o','--outpath',help='Directory for the output, default is the current working directory',default='./')
args = parser.parse_args()
if args.path.endswith('/'):
    pathway=args.path[:-1]
else:
    pathway=args.path
if not args.outpath.endswith('/'):
    outpathway=args.outpath+'/'
else:
    outpathway=args.outpath
if not os.path.isdir(pathway):
    raise Exception('The input pathway does not exist or is not a directory')
if not os.path.exists(outpathway):
    os.mkdir(outpathway)

grabOutput(pathway)
modifyIntervalFile(args.interval,args.number)
findIntersect(args.interval,args.number,outpathway)

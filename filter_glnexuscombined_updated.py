import os
from argparse import ArgumentParser
from re import I


def checkParents(parent,gq_filter,depth_filter,dp_i,gq_i,ad_i):
    #parent=parent.split(',')
    #print(parent)
    depth=int(parent[dp_i])
    allref=int(parent[ad_i].split(',')[1])
    gq=int(parent[gq_i])

    if depth >= depth_filter and allref ==0 and gq>gq_filter:
    #    print(depth,allref,gq,'parent')
        return(1)
    else:
        return(0)

def checkKid(kid,gq_filter,depth_filter,dp_i,gq_i,ad_i):
    depth=int(kid[dp_i])
    ref=int(kid[ad_i].split(',')[1])
    #total=int(depth+ref)
    gq =int(kid[gq_i])

    if depth>=depth_filter and gq > gq_filter:
        #total=int(depth)+int(ref)
        checkme=float(ref)/float(depth)
        #print(checkme)
        if float(checkme) > .25:
        #    print(depth,gq,checkme,'kid')
            return(1)
        else:
            return(0)
    else:
        return(0)

def checkPos(checkme,first,second,third):
    if checkme==first:
        return(0)
    elif checkme==second:
        return(1)
    elif checkme==third:
        return(2)
#flie= sys.argv[1]
#file='HG00405.glnexus_denovo_only_combined_intersection.vcf'
parser = ArgumentParser(description='Grab info for filter')
parser.add_argument('file',help='File to run')
parser.add_argument('f',help='Father')
parser.add_argument('m',help='Mother')
parser.add_argument('c',help='Child')
parser.add_argument('first',help='P1 from combined')
parser.add_argument('second',help='P2 from combined')
parser.add_argument('third',help='P3 from combined')
parser.add_argument('gq_filter',help='GQ filter value',type=int)
parser.add_argument('depth_filter',help='Depth filter value',type=int)
args = parser.parse_args()
holdme=[]
header=[]
file=args.file
name=file.strip().split('.glnexus.family.combined_intersection.vcf')[0]
fi=checkPos(args.f,args.first,args.second,args.third)
mi=checkPos(args.m,args.first,args.second,args.third)
ci=checkPos(args.c,args.first,args.second,args.third)


with open(file) as input:

    for line in input:

        if not line.startswith('#'):
            data=line.strip().split('\t')
            samples=data[9:]
            father=samples[fi].split(':')
            mother=samples[mi].split(':')
            kid=samples[ci].split(':')
            index=data[8].split(':')
            dp=0
            ad=0
            gq=0
            #GT:DP:AD:GQ:PL:RNC
            i=0
            while i < len(index):
                if index[i] =='DP':
                    dp=i
                elif index[i]=='AD':
                    ad=i
                elif index[i]=='GQ':
                    gq=i
                i+=1
            #print(dp,ad,gq)
            m=0
            f=0
            #print(father,data)
           # print(father[ad])
            if len(father[ad].split(','))==2 and len(mother[ad].split(','))==2 and len(kid[ad].split(','))==2:
                #print('hello')
                f=checkParents(father,int(args.gq_filter),int(args.depth_filter),dp,gq,ad)
                m=checkParents(mother,int(args.gq_filter),int(args.depth_filter),dp,gq,ad)
                k=checkKid(kid,int(args.gq_filter),int(args.depth_filter),dp,gq,ad)
                #print(f,m)
            if f==1 and m==1 and k==1:
                holdme.append(line)
        else:
            header.append(line)
with open(name+'.glnexus.family.combined_intersection_filtered_gq_'+str(args.gq_filter)+'_depth_'+str(args.depth_filter)+'.vcf','w') as out:
    for thing in header:
        out.write(thing)
    for thing in holdme:
        out.write(thing)


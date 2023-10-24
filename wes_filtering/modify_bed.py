from argparse import ArgumentParser

parser = ArgumentParser(description='Grab info for filter')

parser.add_argument('-i','--interval',help='path to interval file')
parser.add_argument('-c','--change',help='number to change cut bedfile by, defaults to 40',default=40,type=int)
args = parser.parse_args()

with open("mod_bedfile/changed_bedfile_pm"+str(args.change)+'.bed','w') as out:
    with open(args.interval) as input:
        for line in input:
            data=line.strip().split('\t')
            if int(data[2])-args.change > int(data[1])+args.change:
                out.write(data[0]+'\t'+str(int(data[1])+args.change)+'\t'+str(int(data[2])-args.change)+'\n')
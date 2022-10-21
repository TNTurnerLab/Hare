#!/bin/python
#Naive inheritance in trios
#Tychele N. Turner, Ph.D.
#December 19, 2018

import sys
import gzip

try:
       fname = sys.argv[1]
except:
       sys.stderr.write("Usage: \n")
       sys.exit(1)

class VariantReader(object):
        def __init__(self, filename):
                self.fh = gzip.GzipFile(filename, 'r')
        def __iter__(self):
                return self
        def next(self):
                while True:
                        line = self.fh.readline()
                        if line.startswith( '#' ):
                                continue
                        if line == "":
                                self.fh.close()
                                raise StopIteration
                        line = line[:-1]
                        return Variant(line.split('\t'))

class Variant(object):
        def __init__(self, row):
                self.chrom = row[0]
                self.end = row[1]
                self.id = row[2]
                self.ref = row[3]
                self.alt = row[4]
                self.qual = row[5]
                self.filter = row[6]
                self.info = row[7]
                self.format = row[8]
                self.father = row[9]
                self.mother = row[10]
                self.proband = row[11]
#                self.sibling = row[12]
#                if not self.chrom.startswith("chr"):
#                        self.chrom = 'chr' + self.chrom
                if (self.father.split(':')[0] == '0/0' or self.father.split(':')[0] == '0|0') and (self.mother.split(':')[0] == '0/0' or self.mother.split(':')[0] == '0|0') and (self.proband.split(':')[0] == '0/1' or self.proband.split(':')[0] == '1/1' or self.proband.split(':')[0] == '1/0' or self.proband.split(':')[0] == '0|1' or self.proband.split(':')[0] == '1|1' or self.proband.split(':')[0] == '1|0'):
                        self.info = self.info + ';TRANSMITTED=no;INH=denovo_pro'
                else:
                        self.info = self.info + ';TRANSMITTED=unknown;INH=unknown'
        def __repr__(self):
                return "\t".join([self.chrom, str(self.end), self.id, self.ref, self.alt, self.qual, self.filter, self.info, self.format, self.father, self.mother, self.proband])

vcfFile = gzip.GzipFile(fname, "r")
for line in vcfFile:
	if line.startswith('##FORMAT'):
		continue
	elif line.startswith('##'):
		print line,

print '##INFO=<ID=TRANSMITTED,Number=A,Type=String,Description="This tells whether or not the allele was transmitted from the parents">'
print '##INFO=<ID=INH,Number=A,Type=String,Description="inheritance pattern assuming the rules of autosomal inheritance">'

vcfFile = gzip.GzipFile(fname, "r")
for line in vcfFile:
        if line.startswith('##FORMAT'):
                print line,

vcfFile = gzip.GzipFile(fname, "r")
for line in vcfFile:
	if line.startswith('#C'):
		print line,

for variant in VariantReader(fname):
        print variant

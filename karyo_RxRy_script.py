"""
Version: 	
Author:		Kyriaki Anastasiadou; Pontus Skoglund
Contact: 	kyriaki.anastasiadou@crick.ac.uk; pontus.skoglund@gmail.com
Date: 		12/01/2023
Citation: 	TBD
Usage:		BAM input (suggested -q 30), estimation of sequencing reads aligning to single chromosomes (chrX, chrY, chr21) over autosomal baseline, sex inference.
Example:	Output:		[Total number of alignments in input] [Rx] [SE for Rx] [Ry] [SE for Ry] [Inferred sex]
"""

import sys
import math

from optparse import OptionParser


usage = "usage: %prog [options] <SAM formatted data from stdin>"
parser = OptionParser(usage=usage)
parser.add_option("--chrXname", action="store", type="string", dest="chrXname",help="Identifier for the X chromosome in the SAM input (use if different than chrX, X etc)",default="X")
parser.add_option("--chrYname", action="store", type="string", dest="chrYname",help="Identifier for the Y chromosome in the SAM input (use if different than chrY, Y etc)",default="Y")
parser.add_option("--chr1name", action="store", type="string", dest="chr1name",help="Identifier for chromosome 1 in the SAM input (use if different than chr1, 1 etc)",default="1")
parser.add_option("--chr2name", action="store", type="string", dest="chr2name",help="Identifier for chromosome 2 in the SAM input (use if different than chr2, 2 etc)",default="2")
parser.add_option("--chr3name", action="store", type="string", dest="chr3name",help="Identifier for chromosome 3 in the SAM input (use if different than chr3, 3 etc)",default="3")
parser.add_option("--chr4name", action="store", type="string", dest="chr4name",help="Identifier for chromosome 4 in the SAM input (use if different than chr4, 4 etc)",default="4")
parser.add_option("--chr5name", action="store", type="string", dest="chr5name",help="Identifier for chromosome 5 in the SAM input (use if different than chr5, 5 etc)",default="5")
parser.add_option("--chr6name", action="store", type="string", dest="chr6name",help="Identifier for chromosome 6 in the SAM input (use if different than chr6, 6 etc)",default="6")
parser.add_option("--chr7name", action="store", type="string", dest="chr7name",help="Identifier for chromosome 7 in the SAM input (use if different than chr7, 7 etc)",default="7")
parser.add_option("--chr8name", action="store", type="string", dest="chr8name",help="Identifier for chromosome 8 in the SAM input (use if different than chr8, 8 etc)",default="8")
parser.add_option("--chr9name", action="store", type="string", dest="chr9name",help="Identifier for chromosome 9 in the SAM input (use if different than chr9, 9 etc)",default="9")
parser.add_option("--chr10name", action="store", type="string", dest="chr10name",help="Identifier for chromosome 10 in the SAM input (use if different than chr10, 10 etc)",default="10")
parser.add_option("--chr11name", action="store", type="string", dest="chr11name",help="Identifier for chromosome 11 in the SAM input (use if different than chr11, 11 etc)",default="11")
parser.add_option("--chr12name", action="store", type="string", dest="chr12name",help="Identifier for chromosome 12 in the SAM input (use if different than chr12, 12 etc)",default="12")
parser.add_option("--chr13name", action="store", type="string", dest="chr13name",help="Identifier for chromosome 13 in the SAM input (use if different than chr13, 13 etc)",default="13")
parser.add_option("--chr14name", action="store", type="string", dest="chr14name",help="Identifier for chromosome 14 in the SAM input (use if different than chr14, 14 etc)",default="14")
parser.add_option("--chr15name", action="store", type="string", dest="chr15name",help="Identifier for chromosome 15 in the SAM input (use if different than chr15, 15 etc)",default="15")
parser.add_option("--chr16name", action="store", type="string", dest="chr16name",help="Identifier for chromosome 16 in the SAM input (use if different than chr16, 16 etc)",default="16")
parser.add_option("--chr17name", action="store", type="string", dest="chr17name",help="Identifier for chromosome 17 in the SAM input (use if different than chr17, 17 etc)",default="17")
parser.add_option("--chr18name", action="store", type="string", dest="chr18name",help="Identifier for chromosome 18 in the SAM input (use if different than chr18, 18 etc)",default="18")
parser.add_option("--chr19name", action="store", type="string", dest="chr19name",help="Identifier for chromosome 19 in the SAM input (use if different than chr19, 19 etc)",default="19")
parser.add_option("--chr20name", action="store", type="string", dest="chr20name",help="Identifier for chromosome 20 in the SAM input (use if different than chr20, 20 etc)",default="20")
parser.add_option("--chr21name", action="store", type="string", dest="chr21name",help="Identifier for chromosome 21 in the SAM input (use if different than chr21, 21 etc)",default="21")
parser.add_option("--chr22name", action="store", type="string", dest="chr22name",help="Identifier for chromosome 22 in the SAM input (use if different than chr22, 22 etc)",default="22")
parser.add_option("--digits", action="store", type="int", dest="digits",help="Number of decimal digits in R_y output",default=4)
parser.add_option("--noheader", action="store_true", dest="noheader",help="Do not print header describing the columns in the output",default=False)
parser.add_option("--idxstats", action="store_true", dest="idxstats",help="Input is from samtools idxstats",default=False)
(options, args) = parser.parse_args()


def binomialSE(estimate,totalnumber):
        return math.sqrt((estimate*(1.0-estimate))/totalnumber)

chrYcount=0
chrXcount=0
chr1count=0
chr2count=0
chr3count=0
chr4count=0
chr5count=0
chr6count=0
chr7count=0
chr8count=0
chr9count=0
chr10count=0
chr11count=0
chr12count=0
chr13count=0
chr14count=0
chr15count=0
chr16count=0
chr17count=0
chr18count=0
chr19count=0
chr20count=0
chr21count=0
chr22count=0
totalcount=0


for line in sys.stdin:
	if line[0] == '@':continue
	col=line.split()
	totalcount += 1
	chromosome=col[2].lstrip('chr')
	

	
	if options.chrYname == chromosome: chrYcount += 1
	elif options.chrXname == chromosome: chrXcount += 1
	elif options.chr1name == chromosome: chr1count += 1
	elif options.chr2name == chromosome: chr2count += 1
	elif options.chr3name == chromosome: chr3count += 1
	elif options.chr4name == chromosome: chr4count += 1
	elif options.chr5name == chromosome: chr5count += 1
	elif options.chr6name == chromosome: chr6count += 1
	elif options.chr7name == chromosome: chr7count += 1
	elif options.chr8name == chromosome: chr8count += 1
	elif options.chr9name == chromosome: chr9count += 1
	elif options.chr10name == chromosome: chr10count += 1
	elif options.chr11name == chromosome: chr11count += 1
	elif options.chr12name == chromosome: chr12count += 1
	elif options.chr13name == chromosome: chr13count += 1
	elif options.chr14name == chromosome: chr14count += 1
	elif options.chr15name == chromosome: chr15count += 1
	elif options.chr16name == chromosome: chr16count += 1
	elif options.chr17name == chromosome: chr17count += 1
	elif options.chr18name == chromosome: chr18count += 1
	elif options.chr19name == chromosome: chr19count += 1
	elif options.chr20name == chromosome: chr20count += 1
	elif options.chr21name == chromosome: chr21count += 1
	elif options.chr22name == chromosome: chr22count += 1


#Chr1 to Chr22 sum, excluding 13, 18 and 21

total_chr=chr1count+chr2count+chr3count+chr4count+chr5count+chr6count+chr7count+chr8count+chr9count+chr10count+chr11count+chr12count+chr14count+chr15count+chr16count+chr17count+chr19count+chr20count+chr22count

#compute R_y

Ry=1.0*chrYcount/(total_chr)
ySE=binomialSE(Ry,(total_chr))

#compute R_x

Rx=1.0*chrXcount/(total_chr)
xSE=binomialSE(Rx,(total_chr))

#compute R_21

R21=1.0*chr21count/(total_chr)
SE_21=binomialSE(R21,(total_chr))

#Define values

XYmeanRy = 0.00246572
XYmeanRx = 0.02969922
XXmeanRy = 0.0000318
XXmeanRx = 0.05877203

XYlowerRy = 0.0019
XYupperRy = 0.0032
XYlowerRx = 0.020
XYupperRx = 0.037
XXlowerRy = 0
XXupperRy = 0.0005
XXlowerRx = 0.047
XXupperRx = 0.063

#Calculate Z score - optional

zXX = (Rx - XXmeanRx)/xSE
zXY = (Ry - XYmeanRy)/ySE

#Assign sex

sex='NA'
if chr1count <= 100:
	sex = 'Not_Assigned_low_coverage'
else:
	if XYlowerRy <= Ry <= XYupperRy:
		if XYlowerRx <= Rx <= XYupperRx:
			sex = 'XY'
		elif XXlowerRx <= Rx <= XXupperRx:
			sex = 'XXY'
		elif Rx <= XYlowerRx:
			sex = 'Consistent_with_XY'
	elif XXlowerRy <= Ry <= XXupperRy:
		if XXlowerRx <= Rx <= XXupperRx:
			sex = 'XX'
		elif XYlowerRx <= Rx <= XYupperRx:
			sex = 'XO'
		elif Rx <= XXlowerRx or XXupperRx <= Rx:
			sex = 'Consistent_with_XX'
	elif XYlowerRy <= (Ry-ySE) <= XYupperRy or XYlowerRy <= (Ry+ySE) <= XYupperRy:
		if XYlowerRx <= (Rx-xSE) <= XYupperRx or XYlowerRx <= (Rx+xSE) <= XYupperRx:
			sex = 'Consistent_with_XY'
		else:
			Sex = 'Contamination'
	elif XXlowerRy <= (Ry-ySE) <= XXupperRy or XXlowerRy <= (Ry+ySE) <= XXupperRy:
		if XXlowerRx <= (Rx-xSE) <= XXupperRx or XXlowerRx <= (Rx+xSE) <= XXupperRx:
			sex = 'Consistent_with_XX'
		else: 
			sex = 'Contamination' 

	else:
     		sex = 'Contamination'



print '\t'.join([str(total_chr),str(Rx),str(xSE),str(Ry),str(ySE),str(sex)])

exit(0)

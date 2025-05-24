import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
import matplotlib.patches as mplpatches
from operator import itemgetter
import argparse
import pysam

"""
readCoverage.py takes a BAM file as input and outputs a visual of where reads map to (.png).

Example:
python3 readCoverage.py -b <BAM file> -o <output.bed>
"""


parser = argparse.ArgumentParser()
parser.add_argument('--bam','-b',default="None",type=str,action='store',help='Path to BAM file goes here')
parser.add_argument('--out','-o',default="None",type=str,action='store',help='Path to PNG output goes here')

bamInput=parser.parse_args().bam
output=parser.parse_args().out

# convert BAM file into a Dataframe with read ID, SAMflag, mapped chr name, CIGAR string,
# alignment score, NM score (edit distance)
save = pysam.set_verbosity(0)
bam_data = []
alignments = pysam.AlignmentFile(bamInput,"rb")
pysam.set_verbosity(save)

for read in alignments:
    readName = str(read.query_name) # read id
    derivedChr = readName.split(':')[0] # chromosome the read was derived from
    startRegion = int(readName.split(':')[1].split('-')[0])
    endRegion = int(readName.split(':')[1].split('-')[1].split('_')[0])
    readStart = int(readName.split('_')[2])
    derivedStart = readStart + startRegion
    derivedEnd = readStart + endRegion
    cigarString = read.cigarstring
    flag = read.flag # SAM flag
    cigarSize = read.infer_query_length() # read size inferred from cigar string
    mappedChr = alignments.get_reference_name(read.reference_id) # mapped chromosome name
    mappedStart = int(read.reference_start) + 2 # mapped coordinate
    mappedEnd = int(read.reference_end) if read.reference_end is not None else 'N/A'
    mapQ = int(read.mapping_quality)
    nmScore = read.get_tag("NM") if read.has_tag("NM") else 'N/A'
    asScore = int(read.get_tag("AS")) if read.has_tag("AS") else 'N/A'
    bam_data.append([readName,derivedChr,readStart,mappedChr,flag,mapQ,cigarSize,cigarString,derivedStart,derivedEnd,mappedStart,mappedEnd,nmScore,asScore])

colNames = ['readID','derivedChr','readStart','mappedChr','flag','MapQ','cigarSize','cigarString','derivedStart','derivedEnd','mappedStart','mappedEnd','nmScore','asScore']

df = pd.DataFrame(bam_data,columns=colNames)
# chr = df.loc[0,mappedChr]
readSize=df.loc[0,'cigarSize']
derivedStartCoord = df.loc[0, 'derivedStart']
derivedEndCoord = df.loc[0, 'derivedEnd']

# obtain uniquely mapped reads
value_counts = df['readID'].value_counts()
unique_mapped_values = value_counts[value_counts == 1].index.tolist()
unique_mapped_df = df[df['readID'].isin(unique_mapped_values)& (df['flag'] != '4')]
unique_mapped_df['mappedStart'] = unique_mapped_df['mappedStart'].astype(int)

# mapped coordinates of uniquely mapped reads
uniqueList = list(zip(unique_mapped_df['mappedStart'], unique_mapped_df['cigarSize']))
uniqueList = sorted(uniqueList, key=lambda x: x[1])

# MULTIMAPPED READS
# plots just the multimapped reads with three or less top alignments
# Extract values that appear more than once
multi_mapped_values = value_counts[value_counts > 1].index.tolist()
# Filter the DataFrame to get the multi-mapped entries
multi_mapped_df = df[df['readID'].isin(multi_mapped_values)& (df['flag'] != '4')]
# filtering for just the TOP multi-mapped reads
topAlignmentsDF = multi_mapped_df[(multi_mapped_df['asScore'] == 2*readSize) & ( multi_mapped_df['nmScore'] == 0)]
topReadCounts = topAlignmentsDF['readID'].value_counts()
# Get a list of alignments that appear 2 times or less
reads_to_keep = topReadCounts[topReadCounts <= 2].index
# Use the isin function to filter the DataFrame for top alignments
top2AlignmentsDF= topAlignmentsDF[topAlignmentsDF['readID'].isin(reads_to_keep)]
top2AlignmentsDF['mappedStart'] = top2AlignmentsDF['mappedStart'].astype(int)
# mapped coordinates of uniquely mapped reads
top2AlignsList = list(zip(top2AlignmentsDF['mappedStart'], top2AlignmentsDF['cigarSize']))
top2AlignsList = sorted(top2AlignsList, key=lambda x: x[1])

Purple= '#8B79A5'
Green ='#95BCA5'
darkGreen ='#325b38'
Blue ='#89B0D0'
Pink='#FFB6C1'

finalList = []

for read in uniqueList:
    read_start = read[0]
    read_end = read_start + read[1]
    block_start = read_start
    block_width = read[1]
    finalList.append([read_start, read_end, block_start, block_width,Purple,False])

for read in top2AlignsList:
    read_start = read[0]
    read_end = read_start + read[1]
    block_start = read_start
    block_width = read[1]
    finalList.append([read_start, read_end, block_start, block_width,darkGreen,False])

figureWidth=8
figureHeight=5
plt.figure(figsize=(figureWidth,figureHeight))
panelHeight=4
panelWidth=8
panel2 = plt.axes([0.1/figureWidth,1.7/figureHeight,panelWidth/figureWidth,panelHeight/figureHeight])

iBlue=(88/255,85/255,120/255)

finalList.sort(key=itemgetter(0))
countDict = {}
count=0

for ypos in range(1,len(finalList)):
    lastVal = derivedStartCoord  # HARDCODE
    for read in finalList:
        gstart,gend,blockstarts,blockwidths,readType,plotted=read[0],read[1],read[2],read[3],read[4],read[5]
        if plotted != True: # if not plotted yet, check if start coord is greater than curr end
            if gstart > lastVal: # if greater -> plot read
                # print(f'plotting rectangles {count}')
                rectangle = mplpatches.Rectangle((blockstarts, ypos),
                                                blockwidths,0.5,
                                                facecolor=readType,
                                                edgecolor='black',
                                                linewidth=0.05)
                panel2.add_patch(rectangle)
                lastVal = gend
                read[5] = True
                count+=1

panel2.set_xlim(derivedStartCoord - 5000,derivedEndCoord + 5000)  # HARDCODE
panel2.set_ylim(-1,350)
panel2.set_xlabel('Genomic Coordinate (Mb)')
# panel2.set_xticks(['57','58','59','60','61','62'])
# plt.title(f'{i} Read Coverage for {j}bp read size')
plt.savefig(output, dpi=2400)
plt.clf()
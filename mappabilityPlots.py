import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
import matplotlib.patches as mplpatches
from operator import itemgetter
import pysam
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--bam','-b',default="None",type=str,action='store',help='Path to BAM folder goes here')
parser.add_argument('--chr','-c',default="None",type=str,action='store',help='chr name goes here')
# parser.add_argument('--start','-s',default=0,type=int,action='store',help='region start coordinate goes here')
parser.add_argument('--outputFile','-o',default='None',type=str,action='store',help='output file goes here')

bamFolder=parser.parse_args().bam
chrName=parser.parse_args().chr
# startCoord=parser.parse_args().start

readLengths = [
            1000,
            5000,
            10000,
            15000,
            20000,
            40000,
            60000,
            80000,
            100000,
            200000,
            300000,
            ]

multi_1,multi_2,multi_3 = [],[],[]
uniqueList,multiMappedSameChr,multiMappedDiffChr,unmappedList,completeAlignmentList = [],[],[],[],[]
for j in readLengths:
    print(f'Analyzing {j} bp reads for {chrName}')

    # parse BAM file and add read ID, samFlag, CIGAR, AS, NM to dataframe
    bamPath = f"{bamFolder}/{j}.bam"
    bam_data = []
    save = pysam.set_verbosity(0)
    bamfile = pysam.AlignmentFile(bamPath,"rb")
    pysam.set_verbosity(save)
    for read in bamfile:
        qname = read.query_name
        flag = read.flag
        rname = bamfile.get_reference_name(read.reference_id)
        cigar = read.cigarstring
        as_score = read.get_tag("AS") if read.has_tag("AS") else None
        nm_score = read.get_tag("NM") if read.has_tag("NM") else None
        bam_data.append([qname,flag,rname,cigar,as_score,nm_score])
    df = pd.DataFrame(bam_data)

    # UNMAPPED READS
    # defined by SAMflag=4
    unmapped_reads = df[df[1] == 4]
    unmappedCount = unmapped_reads.shape[0]

    # drop unmapped reads
    df = df[df[1] != 4]

    # Count the occurrences of QNAMES in the first column (read IDs)
    value_counts = df[0].value_counts()

    # UNIQUELY MAPPED READS
    # First, extract read IDs that appear just once
    unique_mapped_values = value_counts[value_counts == 1].index.tolist()
    # Filter all reads to get just uniquely mapped reads dataframe
    unique_mapped_df = df[df[0].isin(unique_mapped_values)]
    # split read ID to retrieve the derived coordinate of the read
    # splitcolumns = unique_mapped_df[0].str.split('_', expand=True)
    # splitcolumns.rename(columns={2: 'derivedCoord'},inplace=True)
    # # create new column with the derived coordinates of each read to unique reads dataframe
    # unique_mapped_df = pd.concat([unique_mapped_df, splitcolumns['derivedCoord']],axis=1)
    # unique_mapped_df['derivedCoord'] = unique_mapped_df['derivedCoord'].astype(int)
    # unique_mapped_df['derivedCoord'] = unique_mapped_df['derivedCoord'] + startCoord - 1
    # # uniquely mapped reads defined by (mapped coordinate == derived coordinate) and maps to the correct haplotype
    # unique_mapped_df = unique_mapped_df[(unique_mapped_df[3] == unique_mapped_df['derivedCoord']) & (unique_mapped_df[2]==chrName)]
    # # num of uniquely mapped reads
    uniquesCount = unique_mapped_df.shape[0]

    # MULTIMAPPED READS
    # Extract read IDs that appear more than once
    multi_mapped_values = value_counts[value_counts > 1].index.tolist()
    # Filter all reads to get the multi-mapped reads
    multi_mapped_df = df[df[0].isin(multi_mapped_values)]
    # is this necessary (yes it is)
    # reads that map to the correct chromosome
    sameChrDF = multi_mapped_df[multi_mapped_df[2].str.startswith(chrName[:5])]
    # reads that map to a different chromosome
    diffChrDF = multi_mapped_df[~multi_mapped_df[2].str.startswith(chrName[:5])] 
    
    # exclude reads that are in the diffChr dataframe
    sameChrDF = sameChrDF[~sameChrDF[0].isin(diffChrDF[0])]
    # filtering for just the top multi aligned reads
    completeAlignment = sameChrDF[(sameChrDF[3] == f'{j}M') & ( sameChrDF[4] == 2*j) & (sameChrDF[5] == 0)]

    read_counts = completeAlignment[0].value_counts()

    # Get a list of reads that appear 3 times or less
    readsTopAligns = read_counts[read_counts <= 2].index
    reads3Aligns = read_counts[read_counts == 3].index
    reads2Aligns = read_counts[read_counts == 2].index
    reads1Aligns = read_counts[read_counts == 1].index

    # Use the isin function to filter the DataFrame for top alignments
    readsTopAlignsDF= completeAlignment[completeAlignment[0].isin(readsTopAligns)]
    reads3AlignDF = completeAlignment[completeAlignment[0].isin(reads3Aligns)]
    reads2AlignDF = completeAlignment[completeAlignment[0].isin(reads2Aligns)]
    reads1AlignDF = completeAlignment[completeAlignment[0].isin(reads1Aligns)]

    # filter for just top 2 alignments to same hap
    newDF = reads2AlignDF[reads2AlignDF[2]==chrName]
    newDFList = newDF[0].value_counts()
    reads2AlignsSameHap = newDFList[newDFList == 2].index
    reads2AlignsSameHapDF = reads2AlignDF[reads2AlignDF[0].isin(reads2AlignsSameHap)]
    reads2AlignsSameHapDF[0].nunique()

    # multi_1.append(len(reads1Aligns))
    # multi_2.append(len(reads2Aligns))
    # multi_3.append(len(reads3Aligns))

    multi_1.append(reads1AlignDF[reads1AlignDF[2]==chrName][0].nunique())
    # multi_2.append(len(reads2Aligns))
    multi_2.append(reads2AlignsSameHapDF[0].nunique())

    uniqueList.append(uniquesCount)
    multiMappedSameChr.append(multi_mapped_df[0].nunique() - diffChrDF[0].nunique()-reads1AlignDF[reads1AlignDF[2]==chrName][0].nunique()-reads2AlignsSameHapDF[0].nunique())
    multiMappedDiffChr.append(diffChrDF[0].nunique())
    unmappedList.append(unmappedCount)
    completeAlignmentList.append(readsTopAlignsDF[0].nunique())

# plot bars in stack manner
x = ["1000","5000","10000","15000","20000","40000","60000","80000","100000",
     "200000","300000"
     ]

y1 = np.array(uniqueList) 
y2 = np.array(multi_1) # reads with 1 top alignment (right hap)
y3 = np.array(multi_2) # reads with 2 top alignments (both right hap)
y4 = np.array(multiMappedSameChr) # everything else (3+)
y5 = np.array(multiMappedDiffChr)
y6 = np.array(unmappedList)
# y1 = np.array(uniqueList) 
# y2 = np.array(completeAlignmentList)
# y3 = np.array(multiMappedSameChr)
# y4 = np.array(multiMappedDiffChr)
# y5 = np.array(unmappedList)

LightGreen= "#90EE90"
LimeGreen="#32CD32"
ForestGreen="#228B22"
Purple= '#8B79A5'
Blue ='#89B0D0'
Pink='#FFB6C1'
print("line155")
plt.bar(x, y1, color=Purple)
plt.bar(x, y2, bottom=y1, color=ForestGreen)
plt.bar(x, y3, bottom=y1+y2, color=LimeGreen)
plt.bar(x, y4, bottom=y1+y2+y3, color=LightGreen)
plt.bar(x, y5, bottom=y1+y2+y3+y4, color=Blue)
plt.bar(x, y6, bottom=y1+y2+y3+y4+y5, color=Pink)
plt.xlabel("Read Length [bp]")
plt.ylabel("Read Count")
plt.xticks(rotation=45)

plt.legend(['uniquely mapped reads', 'reads with 1 top alignment (correct hap)', 'reads with 2 top alignments (both same hap)', 'reads with 3+ alignments or 2 top to both haps','reads mapped to diff chr','unmapped'])
plt.title(f"Read Breakdown for {chrName}")

outFile = parser.parse_args().outputFile
plt.savefig(outFile, dpi=300)
plt.clf()

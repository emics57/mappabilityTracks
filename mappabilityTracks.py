import pandas as pd
import argparse
import pysam

"""
mappability.py takes a BAM file as input and outputs the mappable regions as a BED file.

Example:
python3 mappabilityBEDtools.py -b <BAM file> -o <output.bed>
"""

parser = argparse.ArgumentParser()
parser.add_argument('--bam','-b',default="None",type=str,action='store',help='Path to BAM file goes here')
parser.add_argument('--out','-o',default="None",type=str,action='store',help='Path to BED output goes here')

bamInput=parser.parse_args().bam

# convert BAM file into a Dataframe with read ID, SAMflag, mapped chr name, CIGAR string,
# alignment score, NM score (edit distance)
save = pysam.set_verbosity(0)
bam_data = []
alignments = pysam.AlignmentFile(bamInput,"rb")
pysam.set_verbosity(save)
for read in alignments:
    readName = read.query_name # read id
    derivedChr = readName.split(':')[0] # chromosome the read was derived from
    flag = read.flag # SAM flag
    cigarSize = read.infer_query_length() # read size inferred from cigar string
    mappedChr = alignments.get_reference_name(read.reference_id) # mapped chromosome name
    readStart = read.reference_start # mapped coordinate
    readEnd = read.reference_end
    bam_data.append([readName,derivedChr,flag,mappedChr,cigarSize,readStart,readEnd])
dataframe = pd.DataFrame(bam_data)
chr = dataframe.loc[0,1]

def createCoords(df):
    """
    This function takes in a dataframe of read alignments and 
    outputs a list of coordinates that define regions where reads uniquely map to.
    """
    # obtain uniquely mapped reads
    value_counts = df[0].value_counts()
    unique_mapped_values = value_counts[value_counts == 1].index.tolist()
    unique_mapped_df = df[df[0].isin(unique_mapped_values)& (df[2] != 4)]
    unique_mapped_df.loc[:,5] = unique_mapped_df[5].astype(int)
    unique_mapped_df.loc[:,6] = unique_mapped_df[6].astype(int)

    # mapped coordinates of uniquely mapped reads
    mappedCoord = list(zip(unique_mapped_df[5], unique_mapped_df[6]))
    mappedCoord.sort()
    mappedCoord = sorted(mappedCoord, key=lambda x: x[0])

    # define uniquely mapped regions
    tupleList=[]
    start = mappedCoord[0][0]
    end = mappedCoord[0][1]
    numReads=0
    for i in range(len(mappedCoord) - 1):
        # if end coord less than the next largest start coord -> start new region 
        currReadStart = mappedCoord[i][0]
        currReadEnd = mappedCoord[i][1]
        nextReadStart = mappedCoord[i+1][0]
        nextReadEnd = mappedCoord[i+1][0]
        numReads+=1
        if currReadEnd < nextReadStart:
            # set mappable end coordinate to current mapped coordinate + size
            end = currReadEnd # curr read's end coordinate
            tupleList.append((start,end,numReads))
            # set start coordinate to next largest start coordinate
            start = nextReadStart
            numReads=0
    end = mappedCoord[i+1][1]
    numReads+=1
    tupleList.append((start,end,numReads))
    return tupleList

def tuples_to_bed(tuples_list, bed_file_path, chrom):
    """
    This function generates a BED file from a list of tuples of (start coord, end coord)
    """
    with open(bed_file_path, 'w') as bed_file:
        for tup in tuples_list:
            start, end, numReads = tup
            bed_line = f"{chrom}\t{start}\t{end}\n"
            bed_file.write(bed_line)


def main():
    output_bed_file=parser.parse_args().out
    # define mappable regions
    coords = createCoords(dataframe)
    # generate bed file
    tuples_to_bed(coords, output_bed_file, chr)

if __name__ == "__main__":
    main()

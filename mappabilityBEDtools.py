import pandas as pd
import argparse
import pysam
import subprocess

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
size = dataframe.loc[0,4]
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
    unique_mapped_df[5] = unique_mapped_df[5].astype(int)
    unique_mapped_df[6] = unique_mapped_df[6].astype(int)
    
    # mapped coordinates of uniquely mapped reads
    mappedCoord = list(zip(unique_mapped_df[5], unique_mapped_df[6]))
    mappedCoord.sort()
    mappedCoord = sorted(mappedCoord, key=lambda x: x[0])
    return mappedCoord

def tuples_to_bed(tuples_list, bed_file_path, chrom):
    """
    This function generates a BED file from a list of tuples of (start coord, end coord)
    """
    with open(bed_file_path, 'w') as bed_file:
        for tup in tuples_list:
            start, end = tup
            bed_line = f"{chrom}\t{start}\t{end}\n"
            bed_file.write(bed_line)


def main():
    merged_bed_file=parser.parse_args().out
    # define mappable regions
    coords = createCoords(dataframe)
    # generate read Coordinate BED file
    readCoordinatesBED = f"{size}bp.{chr}.readsCoordinates.bed"
    tuples_to_bed(coords, readCoordinatesBED, chr)
    # generate merged BED file
    command = f"bedtools merge -i {readCoordinatesBED} -c 1 -o count > {merged_bed_file}"
    subprocess.run(command, shell=True)

if __name__ == "__main__":
    main()

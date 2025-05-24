import seaborn as sns
from shiny import reactive
from pathlib import Path
import pysam
import pandas as pd
from operator import itemgetter
import matplotlib.patches as mplpatches
import matplotlib.pyplot as plt 
from functools import partial
from shiny.express import input, render, ui
from shiny.ui import page_navbar

ui.page_opts(
    title="Genome Mappability Visualizer",
    page_fn=partial(page_navbar, id="page"), 
    fillable=True
    )

ui.nav_spacer()  # Push the navbar items to the right

footer = ui.input_select(
    "var", "Select variable", choices=["bill_length_mm", "body_mass_g"]
)

@reactive.calc
def readBAMFile():
    # bamFile = Path(__file__).parent / "100000.bam"
    bamFile = input.file1()[0]['datapath']
    print(bamFile)
    alignments = pysam.AlignmentFile(bamFile,"rb")
    bam_data = []
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

    bamDF = pd.DataFrame(bam_data,columns=colNames)

    return bamDF

here = Path(__file__).parent

with ui.nav_panel("About"):
    "Read Categorization for Genome Mappability"
    @render.image
    def backgroundImage():
        img = {"src": here / "imgs/mapCategories.png", "width": "800px"}  
        return img

    # def mappabilityCategorizationImage():


with ui.nav_panel("Upload Data"):
    MAX_SIZE = 50000
    ui.input_file("file1", "Choose a file to upload:", multiple=True)
    ui.input_radio_buttons("type", "Type:", ["BAM", "SAM"])


with ui.nav_panel("Genome Viewer"):
    with ui.navset_card_underline(title="Mappability Coverage Tracks", footer=footer):
        with ui.nav_panel("Plot"):

            @render.plot
            def hist():
                # p = sns.histplot(
                #     df, x=input.var(), facecolor="#007bc2", edgecolor="white"
                # )
                # return p.set(xlabel=None)
                df = readBAMFile()
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
                # fig = plt.figure(figsize=(figureWidth,figureHeight))
                # panelHeight=4
                # panelWidth=8
                # panel2 = fig.add_axes([0.1/figureWidth,1.7/figureHeight,panelWidth/figureWidth,panelHeight/figureHeight])
                fig, panel2 = plt.subplots(figsize=(figureWidth, figureHeight))


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
                panel2.set_ylim(-1,80)
                panel2.set_xlabel('Genomic Coordinate (Mb)')
                return fig

        with ui.nav_panel("Table"):

            @render.data_frame
            def data():
                # return df[["species", "island", input.var()]]
                return readBAMFile()

    

    

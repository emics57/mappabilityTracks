## Background
This project was first conceptualized to understand the mappability of human centromeres at different read sizes, specifically interested in long reads and understanding at what read size (most of) the centromere is mappable. Mappability can be defined as whether unique reads map to a certain region or not. If a region does not have any uniquely mapped reads mapped to it, then it is not considered a mappable region. As read size increases, the number of mappable regions should also increase, and we want to know what is the minimum read size we can be confident in regards to centromere mappability. 

## Methods:
1. Subset the specified regions from the genome assembly (e.g. the centromeres)
2. Simulate centromere reads using [wgsim](https://github.com/lh3/wgsim)
3. Align simulated reads to the entire genome assembly using aligner of choice (minimap2 or winnowmap)
4. Downstream Analysis of the resulting BAM files - categorized which reads were uniquely mapped
5. The regions where the uniquely mapped reads mapped to would be considered 'mappable'. BED files (mappability tracks) of where the mappable regions are are generated

## To obtain mappability tracks of centromeres:
1. Download all the scripts within the ```MappabilityTrackWorkflow folder``` to your HPC server
2. Set up your environment before running the script:

    ```conda create -n cenmap samtools=1.16.1 minimap2=2.26 winnowmap=2.03 pysam=0.21.0 pandas```

    You don't have to activate the conda environment before running the script because the script should activate it for you.
    Be sure to name the conda environment properly as ```cenmap```.

4. Run ```CenMap.sh```

## CenMap.sh Parameters: 
```
-g [required; string] Path to genome assembly (fasta file)
-b [optional; string] Path to BED file of specified regions. Default to generating tracks genomewide, but running genomewide has not been tested, so it is recommended to specify regions.
-n [optional; int] Coverage of simulated reads. Default to generating 30x coverage
-a [optional; string] Type of aligner (Minimap2 or Winnowmap) Default to minimap2. To use winnowmap, specify '-a winnow'
```

## Examples on how to run CenMap.sh:
```
# Uses Minimap2 to map reads at 30x coverage
sbatch CenMap.sh -g hg002v1.0.1.fasta -b hg002v1.0.1.activehor.expanded5MB.bed

# Uses Winnowmap to map reads at 30x coverage
sbatch CenMap.sh -g hg002v1.0.1.fasta -b hg002v1.0.1.activehor.expanded5MB.bed -a winnow
```

## Other Scripts:
```mappabilityTracks.py```: manually finds the mappable regions by iterating through all uniquely mapped read coordinates

```mappabilityBEDtools.py```: outputs uniquely mapped reads coordinates as a BED file and uses [bedtools merge](https://bedtools.readthedocs.io/en/latest/content/tools/merge.html) to merge reads together to get mappable regions

Both of these scripts (should) do the same thing and output the same exact regions in the file. It's just that the latter requires BEDtools. The bedtools version was created as a sanity check for mappabilityTracks.py.

```mappabilityPlots.py```: Outputs a graph of 10,000 simulated, aligned reads categorized by different mappability definitions (uniquely mapped, multi-mapping, unmapped).

```readCoverage.py```: Outputs a visual of where the uniquely mapped and multi-mapped reads map to.

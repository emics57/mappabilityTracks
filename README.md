```mappabilityTracks.py```: manually finds the mappable regions by iterating through all uniquely mapped read coordinates

```mappabilityBEDtools.py```: outputs uniquely mapped reads coordinates as a BED file and uses [bedtools merge](https://bedtools.readthedocs.io/en/latest/content/tools/merge.html) to merge reads together to get mappable regions

Both of these scripts (should) do the same thing and output the same exact regions in the file. It's just that the latter requires BEDtools. The bedtools version was created as a sanity check for mappabilityTracks.py.

### WES sample filtering for high and low confidence regions



Config Variables:
```
{
  "reference": "/path/to/reference",  #Path to reference file
  "famfile": "/path/to/familyfile",  #Path to family file, formatted as FatherID,MotherID,ChildID
  "data_dir": "/path/to/dnv/files", #Path to folder with vcf files
  "out_dir": "/out_dir_name", #Path to outdir
  "cramloc": "/file/of/cram/locations", #File with cram files locations, one per line
  "lower": 0, #Number of allowed lower quality reads found in the parents
  "upper": 0, #Number of allowed high quality reads found in the parents
  "case_specific": "False", #check if low or high quality reads are found
  "interval": "/path/to/capture_region/file",  #capture file used for exons
  "change": 40 #Defines how many basepairs from each end of the capture we will use to set low/high confidence regions
}
```

## Support code for VarDict assessment

### complex

Compare complex indel calling between VarDict, Scalpel and Pindel

Uses a [snakemake](https://bitbucket.org/snakemake/snakemake/wiki/Home) script,
run with:

   cd complex
   snakemake

It will download original VCF calls, process to BED and compare overlaps in the
validation directory.

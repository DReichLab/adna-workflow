This is the software used to run the Reich Lab ancient DNA workflows. 

Workflows are written in the [Workflow Description Language (WDL)](https://software.broadinstitute.org/wdl/) and are run using the [Broad Institute Cromwell](https://github.com/broadinstitute/cromwell) workflow tool. 

Workflows are setup to run on the Harvard Medical School O2 SLURM cluster. Running on other platforms will may require modifications. 

The workflows interact with a Django website and database to track the results of samples across multiple sequencing runs. 

The workflows use numerous external programs:

- [Reich Lab aDNA tools](https://github.com/DReichLab/ADNA-Tools)
- [The Broad Institute Picard Tools](https://broadinstitute.github.io/picard/)
- [samtools](http://www.htslib.org/)
- [htsbox](https://github.com/lh3/htsbox) for MT consensus calling
- contammix contact Philip Johnson <plfj@umd.edu>
- [haplogrep](http://haplogrep.uibk.ac.at/)
- [bwa](https://github.com/lh3/bwa)
- [pmdtools](https://github.com/pontussk/PMDtools)

Workflows are setup to run on the *scratch* (temporary) filesystem, then copy permanent results to the *group* filesystem. 

- demultiplex.wdl - This takes an Illumina sequencer output directory as input and outputs a series of bams named by the index and barcode. This is run once per sequencing run. 

 These bams follow the naming pattern [i5 index]\_[i7 index]\_[p5 barcode]\_[p7 barcode], and are stored on the permanent filesystem. 

 Paired-end reads are merged into single-end reads, requiring some minimum overlap and allowing for some mismatch depending on base quality scores. Adapters are trimmed during this while merging. 

 There are two sets of bams, one aligned to the whole human genome reference hg19, and one aligned to the mitochondrial Reconstructed Sapiens Reference Sequence (RSRS). Bams are filtered to include only reads aligning to the reference. 

- analysis.wdl - This calculates a number of metrics for each bam, both on the Reich Lab set of ~1240k nuclear target data and MT data. It builds bams for each sample based on prior sequencing runs.

- release\_and\_pulldown.wdl - Build release versions of bams with one read group per flowcell lane. Run pulldown to generate pseudo-haploid genotype data. 

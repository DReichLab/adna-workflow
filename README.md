This is the software used to run the Reich Lab ancient DNA workflow. 

This workflow uses the [Broad Institute Cromwell](https://github.com/broadinstitute/cromwell) workflow tool. 

The workflow is setup to run on the Harvard Medical School O2 SLURM cluster. Running on other platforms will may require modifications. 

The workflow interacts with a Django website and database to track the results of samples across multiple sequencing runs. 

This workflow uses numerous external programs:
-[Reich Lab aDNA tools](https://github.com/DReichLab/ADNA-Tools)
-[The Broad Institute Picard Tools](https://broadinstitute.github.io/picard/)
-[samtools](http://www.htslib.org/)
-[htsbox](https://github.com/lh3/htsbox) for MT consensus calling
-contammix contact Philip Johnson <plfj@umd.edu>
-[haplogrep](http://haplogrep.uibk.ac.at/)
-[bwa](https://github.com/lh3/bwa)
-[pmdtools](https://github.com/pontussk/PMDtools)


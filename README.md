# DiffExp
Differential Expression analysis tools for Galaxy bioinformatics platform. Allows analysis of RNASeq and shRNASeq experiments using generalized linear models through Galaxy. Analysis is performed purely through R scripts with xml wrappers for interface to Galaxy.

**The contents of this repository are intended only for review and development purposes, for the functional version of the tools please visit the relevant Galaxy toolshed repository**

### Functioning Tools
#### shrnaseq (hairpinTool)
https://toolshed.g2.bx.psu.edu/view/shians/shrnaseq

shrnaseq tool (marked as hairpintool) performs sequence alignment, feature counting and statistical summary on shRNA or (s)gRNA sequencing experiments. The analysis methods are based on edgeR Bioconductor package.

More details on workflow can be found at http://bioinf.wehi.edu.au/shRNAseq/

#### voom_rnaseq (diffexp)
https://toolshed.g2.bx.psu.edu/view/shians/voom_rnaseq

voom_rnaseq (marked as diffexp) performs statistical summary on tables of feature counts from RNASeq experiments. The analysis methods are based on the limma Bioconductor package with specific use of the voom function.

More details on workflow can be found at http://bioinf.wehi.edu.au/RNAseqCaseStudy/

# scRNA-Analysis

It`s not easy to run Seurat after version 5 released, since many packages such as doubletFinder and matrix have updated their requiremnets.

Here I present the humble way to do scRNA-seq analysis,source code changes have made to doubleFinder(paramSweep_v3, DoubletFinder_v3), and most command were encapsulated in a function.

You could easily adopt loading/filter/de doublets/de batch effect/SingleR by calling the speficic function after slight modifications on url/data names tec.

Most requirements are listed in the annotations.


***General Rule:***


Please run the **cellranger.sh** to decompress the fastq once you have collected necessary files ( which are barcodes.tsv, features.tsv, matrix.mtx.) and then apply **pipeline S4** on them.


A detailed tutorial with manually annotations method would be coming soon.

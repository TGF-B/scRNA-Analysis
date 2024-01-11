**scRNA-Analysis,from FastQ to SingleR**

It`s not easy to run Seurat after version 5 released, since many packages such as doubletFinder and matrix have updated their requirements.

Here I present the humble way to do scRNA-seq analysis,source code changes have made to doubleFinder(paramSweep_v3, DoubletFinder_v3), and most commands were encapsulated in a function.

You could easily adopt loading/filter/de doublets/de batch effect/SingleR by calling the speficic function after slight modifications on url/data names etc.

Most requirements are listed in the annotations.


Please first run the **cellranger.sh** to decompress the fastq , once you have collected necessary files ( which are barcodes.tsv, features.tsv, matrix.mtx.) , apply **pipeline S4** on them.

Jupyter is basically for entertaining(not really

***Updating:***

A detailed tutorial with manually annotations method would be coming soon.

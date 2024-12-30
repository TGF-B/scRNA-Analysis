**scRNA-Analysis,from FastQ to SingleR**

It`s not easy to run Seurat after version 5 released, since many packages such as doubletFinder and matrix have updated their requirements.

Here I present the humble way to do scRNA-seq analysis,source code changes have made to doubleFinder(paramSweep_v3, DoubletFinder_v3), and most commands were encapsulated in a function.

You could easily adopt loading/filter/de doublets/de batch effect/SingleR by calling the speficic function after slight modifications on url/data names etc.

Most requirements are listed in the annotations.


Please first run the **cellranger.sh** to decompress the fastq , once you have collected necessary files ( which are barcodes.tsv, features.tsv, matrix.mtx.) , apply **pipeline S4** on them.

Jupyter is basically for entertaining(not really

***Updating:***

A detailed tutorial with manually annotations method would be coming soon.

A full version updated:from FASTQ to sub-celltype analysis.(04.10.2024) 

![image](https://github.com/TGF-B/scRNA-Analysis/assets/68436923/81cfcb1c-5861-48a0-99e7-41da5953e659)

Figure 1. scRNA-seq of Cd11b+ and Cd11c+ or Cd11b+/Cd11c+ cells revealed that APC cells increased after CDSA.
![image](https://github.com/TGF-B/scRNA-Analysis/assets/68436923/67b7439a-7af4-4bc8-8bb6-085ca4f69324)

Figure 3: pDC(Runx2+) showed to be homing and recirculated subtype pDC.

***Updating:***
A python version from fastq to scVelo was pushed.(12.30.2024)
With automation by GPT and powerful third-party proceeding, single cell analysis from fastq to cell annotation was much easier.
The average time cost was reduced from 2 months to 1 day...
![image](https://github.com/TGF-B/scRNA-Analysis/blob/main/RNA%20volocity.png)
FIgure 3. RNA volocity prediction based on spliced/unspliced genes.


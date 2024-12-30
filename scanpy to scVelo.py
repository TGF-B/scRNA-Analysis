import scanpy as sc
import scvelo as scv
import anndata as ad
import os
import leidenalg
import anndata as ad
import loompy
import numpy as np
import ipywidgets
import gseapy as gp
import numpy as np

# Ensure the directory exists
output_dir = "/media/donaldtangai4s/Yu_Omics2/plots"
os.makedirs(output_dir, exist_ok=True)


# Function to load data from a single sample
def load_sample(sample_path, sample_name):
    adata = sc.read_10x_mtx(
        sample_path,
        var_names='gene_symbols',
        cache=True
    )
    adata.var_names_make_unique()
    adata.obs['sample'] = sample_name
    return adata

# Load all samples
samples = ['CP-1', 'CP-2', 'CP-3', 'Ctrl-1', 'Ctrl-2', 'Ctrl-3']
adatas = []
#0:CP-1,1:CP-2,2:Cp-3,3:Ctrl-1,4:Ctrl-2,5:Ctrl-3

for sample in samples:
    path = f"/home/donaldtangai4s/new/{sample}/outs/filtered_feature_bc_matrix"
    adata = load_sample(path, sample)
    adatas.append(adata)

# Concatenate all samples
adata_combined = ad.concat(adatas, label="sample")
#rename adata_combined.obs.sample accoring to {'CP-1': 0, 'CP-2': 1, 'CP-3': 3,'Ctrl-1': 2, 'Ctrl-2': 4, 'Ctrl-3': 5})
adata_combined.obs['sample'] = adata_combined.obs['sample'].cat.rename_categories({0: 'CP-1', 1: 'CP-2', 2: 'Ctrl-1', 3: 'CP-3', 4: 'Ctrl-2', 5: 'Ctrl-3'})
#adata.obs['sample'].value_counts()
# Calculate QC metrics
sc.pp.calculate_qc_metrics(adata_combined, inplace=True)

# Plot QC metrics
# 'pct_counts_mt' already filtered by cellranger
sc.pl.violin(adata_combined, ['n_genes_by_counts', 'total_counts'],
             jitter=0.4, multi_panel=True)

# Filter cells
sc.pp.filter_cells(adata_combined, min_genes=200)
sc.pp.filter_genes(adata_combined, min_cells=3)

print(f"Number of cells after filtering: {adata_combined.n_obs}")
print(f"Number of genes after filtering: {adata_combined.n_vars}")


# Calculate QC metrics
sc.pp.calculate_qc_metrics(adata_combined, inplace=True)

# Plot QC metrics
sc.pl.violin(adata_combined, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

# Filter cells
sc.pp.filter_cells(adata_combined, min_genes=200)
sc.pp.filter_genes(adata_combined, min_cells=3)


# Normalize the data
sc.pp.normalize_total(adata_combined, target_sum=1e4)
sc.pp.log1p(adata_combined)

# Identify highly variable genes
sc.pp.highly_variable_genes(adata_combined, min_mean=0.0125, max_mean=3, min_disp=0.5)

# Scale the data
sc.pp.scale(adata_combined, max_value=10)
adata=adata_combined

sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True)

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.tsne(adata)

sc.tl.leiden(adata)

sc.pl.umap(adata, color=['leiden', 'sample'],save='umap.pdf')
sc.pl.tsne(adata, color=['leiden', 'sample'],save='tsne.pdf')
#save the umap and tsne fig

sc.tl.rank_genes_groups(adata, 'subtype', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
#save the preprocessed adata AS H5AD file\
adata.write('preprocessed_adata.h5ad')



# Annotate cell types based on marker genes
# This step requires manual inspection and domain knowledge
cell_type_markers = {
    'T cells': ['Cd3d', 'Cd3e'],
    'B cells': ['Cd79a', 'Ms4a1', 'Pecam1'],
    'NK cells': ['Ncr1', 'Eomes', 'Klrb1c','Bcl2l11'],
    'cDC': ['Batf3', 'Ccr7', 'Ccl2'],
    'pDC': ['Siglech', 'Ly6d'],
    'Mono/Mac': ['Csf1r', 'C1qa', 'C1qb'], 
    'Neutrophils': ['Bcl2l11','S100a8'],
    'Mast cells': ['Cpa3'],
    'Epithelial cells': ['Ly6d','Krt8'],
    'Fibroblasts': ['Dcn'],
    'Endothelial cells': ['Pecam1', 'Cldn5'],
}

def annotate_cell_types(adata):
    cell_types = []
    for cluster in adata.obs['leiden'].cat.categories:
        expr = adata[adata.obs['leiden'] == cluster, :].X.mean(axis=0)
        max_score = 0
        cell_type = 'Unknown'
        for ct, markers in cell_type_markers.items():
            score = np.mean([expr[adata.var_names.get_loc(m)] for m in markers if m in adata.var_names])
            if score > max_score:
                max_score = score
                cell_type = ct
        cell_types.append(cell_type)
    cluster_to_cell_type = dict(zip(adata.obs['leiden'].cat.categories, cell_types))
    adata.obs['cell_type'] = adata.obs['leiden'].map(cluster_to_cell_type)


annotate_cell_types(adata)

sc.pl.umap(adata, color='cell_type')

#load the processed_adata.h5ad as adata
adata = ad.read('preprocessed_adata.h5ad')# # Read loom file


# Load spliced and unspliced matrices
def load_loom(sample_name):
    path = f"/media/donaldtangai4s/Yu_Omics2/{sample_name}.loom"
    return ad.read_loom(path)
samples = ['CP-1', 'CP-2', 'CP-3', 'Ctrl-1', 'Ctrl-2', 'Ctrl-3']
loom_adatas = [load_loom(sample) for sample in samples]

# Ensure unique indices for each AnnData object
for i, adata in enumerate(loom_adatas):
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

# Merge loom files
ldata = ad.concat(loom_adatas, label="sample", index_unique="-")
#save the ldata
ldata.write('merged_loom_data.h5ad')

#load the ldata
ldata = ad.read('merged_loom_data.h5ad')
#extract the pDC from adata
pDC_adata = adata[adata.obs['cell_type'] == 'pDC', :].copy()
#erase the mt gene from pDC
pDC_adata = pDC_adata[:, ~pDC_adata.var_names.str.startswith('mt-')]
#find the top markers of pDC by cluster
sc.tl.rank_genes_groups(pDC_adata, 'leiden', method='wilcoxon')
#plot the most different gene of pDC by CP to Ctrl
sc.pl.rank_genes_groups(pDC_adata, n_genes=25, sharey=False, save='rank_genes_groups.pdf')


# Define subtype markers
pDC_subtype_markers = {
    'pDC(Spib+)': ['Spib'],
    'pDC(Wdfy4+)': ['Wdfy4'],
    'pDC(Runx2+)': ['Runx2'],
    'pDC(Irf7+)':['Irf7','Irf8'],
}

# Annotate subtypes
def annotate_subtypes(adata):
    subtypes = []
    for cluster in adata.obs['leiden'].cat.categories:
        cluster_adata = adata[adata.obs['leiden'] == cluster, :]
        expr = cluster_adata.X.mean(axis=0).A1 if isinstance(cluster_adata.X, np.matrix) else cluster_adata.X.mean(axis=0)
        max_score = 0
        subtype = 'Unknown'
        for st, markers in subtype_markers.items():
            markers_present = [m for m in markers if m in adata.var_names]
            if not markers_present:
                continue
            scores = expr[[adata.var_names.get_loc(m) for m in markers_present]]
            score = scores.mean()
            if score > max_score:
                max_score = score
                subtype = st
        subtypes.append(subtype)
    cluster_to_subtype = dict(zip(adata.obs['leiden'].cat.categories, subtypes))
    adata.obs['subtype'] = adata.obs['leiden'].map(cluster_to_subtype)

# Perform Leiden clustering on pDC cells
sc.tl.leiden(pDC_adata)
subtype_markers=pDC_subtype_markers
# Annotate subtypes
annotate_subtypes(pDC_adata)
#save annotated pDC
pDC_adata.write('pDC_annotated.h5ad')


####################### start here #######################

import scanpy as sc
import scvelo as scv
import anndata as ad
import os
import leidenalg
import anndata as ad
import loompy
import numpy as np
import ipywidgets
import gseapy as gp

pDC_adata=ad.read('pDC_annotated.h5ad')

###########################1:Show the celltyp distribution###############

#chenge color: for example
#UMAP with color setting:pDC(Runx2+) as #E8c023,pDC(Spib+)as #A42A80,
#pDC(Wdf4+) as  #1505808, pDC(Irf+) as #A15E52
sc.pl.umap(pDC_adata, color='subtype', 
           palette=['#A15E52','#150580','#E8c023','#A42A80' ],
           save='pDC_subtype.pdf')

#2:load the neceserry spliced info
#create another column in obs as group which separate 0
#split the UMAP by group which is Ctrl Vs CP
sc.pl.umap(pDC_adata, color=('subtype'),save='pDC_group.pdf')
# Merge with processed adata
adata_velocity = scv.utils.merge(pDC_adata, ldata)
#save the metadata
adata_velocity.write('pDC_velocity.h5ad')

#3: calculate the velocity and velocity graph
#load the pDC_velocity.h5ad
adata_velocity = ad.read('pDC_velocity.h5ad')
scv.pp.filter_and_normalize(adata_velocity, min_shared_counts=10, n_top_genes=3000)
scv.pp.moments(adata_velocity, n_pcs=30, n_neighbors=40)

scv.tl.velocity(adata_velocity, mode='stochastic')
scv.tl.velocity_graph(adata_velocity)

################## 4:plot the peudotime analysis#####################################
# Visualize RNA velocity by group Ctrl and CP horizontally
scv.pl.velocity_embedding_stream(adata_velocity, basis='umap', 
                                 color='subtype', palette=['#150580','#E8c023',
                                                           '#A15E52','#A42A80'],
                                                           save='velocity_stream2.pdf')
####################5:Check interested genes distribution################################
#check  Runx2,Spib+,Wdfy4 and Tcf4 in pDC_adata
sc.pl.umap(pDC_adata, color=['Runx2','Spib', 'Wdfy4', 'Ccr7','Irf7'],save='pDC_markers.pdf')

# Perform differential velocity analysis
scv.tl.rank_velocity_genes(adata_velocity, groupby='leiden', min_corr=.3)
df = scv.get_df(adata_velocity, 'rank_velocity_genes/names')
print(df.head())

# Save results
adata.write('combined_analysis_results.h5ad')
adata_velocity.write('velocity_analysis_results.h5ad')

# Perform GSEA pathway analysis
top_genes = pDC_adata.var_names[pDC_adata.var['highly_variable']]
enr = gp.enrichr(gene_list=top_genes.tolist(),
                 gene_sets=['KEGG_2016', 'GO_Biological_Process_2023'],
                 organism='mouse',  # change to 'mouse' if needed
                 outdir='figures')

# Visualize GSEA results
gp.plot.barplot(enr.res2d, title='Enrichment Analysis', ofname='enrich_barplot.pdf')

# %% [markdown]
# SSkin Warts run for data

# %% [markdown]

# %% [markdown]

# %%
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import matplotlib as mpl
import numpy as np
import pandas as pd
import scanpy as sc
import sys
import os
import scanpy_plus as scp

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=80, facecolor='white', color_map='viridis')
sc.logging.print_header()


input_mtx = os.path.join("../../data/", sys.argv[1], "outs", "count",
                         "filtered_feature_bc_matrix")
results_file = 'pbmc10k_scanpy.h5ad'


# %% [markdown]
# read only gex part

# %%
adata = sc.read_10x_mtx(input_mtx, cache = True)
adata = data[:, adata.var['feature_types'] == 'Gene Expression']
adata.describe()

# %%
adata.var_names_make_unique()

#%% [markdown]
# Filter cells
# %%
sc.pp.filter_genes(adata, min_cells=3)
adata

# %% [markdown]

# %%
sc.pl.highest_expr_genes(adata, n_top=20, )


# %%
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
adata.var['rb'] = adata.var_names.str.startswith(('RPS','RPL'))  # annotate the group of ribosomal proteins as 'rb'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pp.calculate_qc_metrics(adata, qc_vars=['rb'], percent_top=None, log1p=False, inplace=True)


# %%
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_rb'], 
             jitter=0.4, multi_panel=True)


# %%
adata.obs.head()

# %% [markdown]
# Detect doublets
# %%
sc.pp.scrublet(adata)


# %%
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')

# %%
sc.pl.scatter(adata, x='total_counts', y='pct_counts_rb')
sc.pl.scatter(adata, x='pct_counts_rb', y='pct_counts_mt')

# %% [markdown]
# There's also a good correlation between the doublet score and number of expressed genes. 

# %%
sc.pl.scatter(adata, x='n_genes_by_counts', y='Doublet_score')

# %% [markdown]
# Intiial filter
# %%
conditions = [
    (adata.obs['Is_doublet'] == True),
    (adata.obs['n_genes_by_counts'] < 500),
    (adata.obs['pct_counts_mt'] > 15),
    (adata.obs['pct_counts_mt'] <= 15) & (adata.obs['n_genes_by_counts'] >=500) & (adata.obs['predicted_doublet'] != True)
]

values = ['Doublet', 'Low_nFeature', 'High_MT', 'Pass']
adata.obs['QC'] = np.select(conditions, values)
adata.obs['QC'] = adata.obs['QC'].astype('category')
adata.obs[['Is_doublet','n_genes_by_counts','pct_counts_mt','QC']].head(10)

# %%
adata.obs['QC'].value_counts()


# %%
sc.pl.violin(adata[adata.obs['QC'] == 'Pass'], ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_rb'], 
             jitter=0.4, multi_panel=True)

# %% [markdown]
# ## PART 2. Normalization and dimensionality reduction.

# %% [markdown]

# %%
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata

# %%
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# %%
adata.var.highly_variable.value_counts()


# %% [markdown]
# Plot the highly variable genes, after and before the normalization.

# %%
sc.pl.highly_variable_genes(adata)

# %%
adata.raw = adata


# %%
adata = adata[:, adata.var.highly_variable]
adata


# %%
#sc.pp.scale(adata, max_value=10)


# %%
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_scatter(adata)

sc.pl.pca_variance_ratio(adata, log=True)


# %%
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)

# %%
sc.tl.umap(adata)

# %%
sc.pl.umap(adata, color= ['n_genes_by_counts', 'total_counts',
                         'pct_counts_mt', 'pct_counts_rb'])

# %%
sc.tl.leiden(adata)

# %%
sc.pl.umap(adata, color=['leiden'])

# %%
s_genes = ['MCM5','PCNA','TYMS','FEN1','MCM7','MCM4','RRM1','UNG','GINS2','MCM6','CDCA7','DTL','PRIM1',
           'UHRF1','CENPU','HELLS','RFC2','POLR1B','NASP','RAD51AP1','GMNN','WDR76','SLBP','CCNE2','UBR7',
           'POLD3','MSH2','ATAD2','RAD51','RRM2','CDC45','CDC6','EXO1','TIPIN','DSCC1','BLM','CASP8AP2',
           'USP1','CLSPN','POLA1','CHAF1B','MRPL36','E2F8']
g2m_genes = ['HMGB2','CDK1','NUSAP1','UBE2C','BIRC5','TPX2','TOP2A','NDC80','CKS2','NUF2','CKS1B',
             'MKI67','TMPO','CENPF','TACC3','PIMREG','SMC4','CCNB2','CKAP2L','CKAP2','AURKB','BUB1',
             'KIF11','ANP32E','TUBB4B','GTSE1','KIF20B','HJURP','CDCA3','JPT1','CDC20','TTK','CDC25C',
             'KIF2C','RANGAP1','NCAPD2','DLGAP5','CDCA2','CDCA8','ECT2','KIF23','HMMR','AURKA','PSRC1',
             'ANLN','LBR','CKAP5','CENPE','CTCF','NEK2','G2E3','GAS2L3','CBX5','CENPA']
cell_cycle_genes = s_genes + g2m_genes
display(len(s_genes))
display(len(g2m_genes))
display(len(cell_cycle_genes))


# %%
adata = adata.raw.to_adata()
adata


# %%
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
len(cell_cycle_genes)

# %%
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)

# %%
sc.pl.umap(adata, color=['Doublet_score','pct_counts_mt','n_genes_by_counts'])

adata = adata[adata.obs['QC'] == 'Pass']
adata

# %%
sc.pl.umap(adata, color=['Doublet_score','pct_counts_mt','n_genes_by_counts'])

# %%
sc.pl.umap(adata, color=['pct_counts_rb','pct_counts_mt','leiden'])


# %%
sc.pl.umap(adata, color=['S_score','G2M_score','leiden'])

# %%
with mpl.rc_context({'figure.figsize': (5.5, 4)}):
    sc.pl.violin(adata, ['S_score','G2M_score'], groupby = 'leiden', stripplot = False)


# %%
sc.pl.umap(adata, color=['n_genes_by_counts','total_counts','leiden'])

# %%
with mpl.rc_context({'figure.figsize': (5.5, 4)}):
    sc.pl.violin(adata, ['n_genes_by_counts','total_counts'], groupby = 'leiden', stripplot = False)

# %% [markdown]
# ### PART 3. Normalization and careful clustering of the cleaned dataset.


# %%
adata


# %%
sc.pp.highly_variable_genes(adata, min_mean = 0.0125, max_mean = 3, min_disp = 0.5)
adata.var.highly_variable.value_counts()

# %% [markdown]
# %%
adata.raw = adata
adata = adata[:, adata.var.highly_variable]
adata

#sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt', 'S_score', 'G2M_score'])

# %% [markdown]
# Calculate PCA and plot the first two components.

# %%
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca(adata, color='CST3')

# %%
sc.pl.pca_variance_ratio(adata, log=True)

# %%
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata)

# %%
sc.tl.leiden(adata,resolution = 0.8)

with mpl.rc_context({'figure.figsize': (8, 8)}):
    sc.pl.umap(adata, color=['leiden'],legend_loc='on data')

# %%
sc.tl.paga(adata)
sc.pl.paga(adata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
sc.tl.umap(adata, init_pos='paga')

# %% [markdown]
# After this, take another look at the updated UMAP. Looks better!

# %%
with mpl.rc_context({'figure.figsize': (8, 8)}):
    sc.pl.umap(adata, color=['leiden'],legend_loc='on data')

# %% [markdown]
# Let's visualize more detailed markers of PBMC cell types:
# 
# | Marker | Cell type | 
# | :-: | :-: |
# | LYZ | monocytes |
# | MS4A1 | B cells |
# | CD8B | CD8 T cells |
# | CCR7 | CD4 naive T cells |
# | IL32 | CD4 memory T cells |
# | NKG7 | natural killer cells |
# | LILRA4 | plasmacytoid dendritic cells |
# | FCER1A | myeloid dendritic cells |
# | PPBP | platelets |

# %%
sc.pl.umap(adata, color=['LYZ','MS4A1','CD8B','CCR7','IL32','NKG7','LILRA4','FCER1A','PPBP'], ncols = 3)

# 
exit(0)

# %% [markdown]
# Let's try some other approaches. For example, we can try k-means clustering with 15, 20, or 25 expected clusters. 

# %%
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
X_pca = adata.obsm['X_pca'] 

kmeans = KMeans(n_clusters=15, random_state=0).fit(X_pca) 
adata.obs['kmeans15'] = kmeans.labels_.astype(str)
kmeans = KMeans(n_clusters=20, random_state=0).fit(X_pca) 
adata.obs['kmeans20'] = kmeans.labels_.astype(str)
kmeans = KMeans(n_clusters=25, random_state=0).fit(X_pca) 
adata.obs['kmeans25'] = kmeans.labels_.astype(str)

sc.pl.umap(adata, color=['kmeans15', 'kmeans20', 'kmeans25'], legend_loc='on data', wspace = 0.25, legend_fontsize=10)

# %% [markdown]
# We can also try hierarchical clustering using Euclidean distance as a metric, and Ward's linkage method. 

# %%
from sklearn.cluster import AgglomerativeClustering

cluster = AgglomerativeClustering(n_clusters=15, affinity='euclidean', linkage='ward')
adata.obs['hclust15'] = cluster.fit_predict(X_pca).astype(str)

cluster = AgglomerativeClustering(n_clusters=20, affinity='euclidean', linkage='ward')
adata.obs['hclust20'] = cluster.fit_predict(X_pca).astype(str)

cluster = AgglomerativeClustering(n_clusters=25, affinity='euclidean', linkage='ward')
adata.obs['hclust25'] = cluster.fit_predict(X_pca).astype(str)

sc.pl.umap(adata, color=['hclust15', 'hclust20', 'hclust25'], legend_loc='on data', wspace = 0.25, legend_fontsize=10)

# %% [markdown]
# In general, k-means and hierarchical clusters seem to poorly represent the communities, easily shifting from under- to over-clustering. At the same time, one of the clusters from our initial Leiden clustering is always split into platelet/DC sub-populations. This indirectly confirms our previous conclusions. Thus, let's manually split this cluster using well-known markers of platelets.

# %%
adata.obs['leiden'].value_counts()

# %% [markdown]
# One of the clusters is clearly divided into two by the following markers (PPBP is a marker of platelets, and FCER1A is a marker of myeloid DCs)

# %%
sc.pl.umap(adata, color=['PPBP','FCER1A','leiden'])

# %% [markdown]
# Let's find which cells of cluster 10 are platelets. The following expression shows a useful way to slice `scanpy` object on multiple gene expression values. The resulting subset should contain 25 cells.

# %%
platelets = adata[(adata[:,'GP1BB'].X > 1) & (adata[:,'PPBP'].X > 1) & (adata[:,'PF4'].X > 1), :]
platelets

# %% [markdown]
# We can get the indexes of putative platelets from cluster 10 as follows:

# %%
platelets.obs_names

# %% [markdown]
# The following commands allow us to add another category to the 'leiden' column of the metadata, and set the selected cells to belong to the new cluster:

# %%
adata.obs['leiden'] = adata.obs['leiden'].cat.add_categories('16')
adata.obs.loc[lambda df: df.index.isin(platelets.obs.index),'leiden'] = '16'

# %% [markdown]
# Finally, let's make sure the clustering looks right: 

# %%
with mpl.rc_context({'figure.figsize': (8, 8)}):
    sc.pl.umap(adata, color=['leiden'],legend_loc='on data')

# %% [markdown]
# ### PART 4. Differential expression and marker selection.

# %% [markdown]
# After we have fixed the clusters, let's make gene markers using Wilcoxon test. Score plotting for each cluster allows to estimate the relative "goodness" of a marker in a quick fashion.

# %%
adata

# %%
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=15, sharey=False)

# %% [markdown]
# While some clusters have uninformative markers (e.g. cluster 0 mostly shows ribosomal proteins), most others identify what seems to be legitimate subpopulations of lymphocytes. Let's define a set of markers that are known from the literature or have been identified in our marker selection process:

# %%
marker_genes = ['IL7R', 'CD79A', 'MS4A1', 'CD8A', 'CD8B', 'LYZ', 'CD14',
                'LGALS3', 'S100A8', 'GNLY', 'NKG7', 'KLRB1', 'IGLC1','IL32',
                'FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'PPBP','LILRA4']

# %% [markdown]
# Let's take a closer look at top 5 ranked markers for each cluster: 

# %%
pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5)

# %% [markdown]
# Seaborn (split violin) plots of gene expression in a given cluster vs the rest of cells can help easily identify more robust cluster markers. In the example below, ribosomal genes should not be used as markers, while CD8A and CD8B represent *bona fide* markers of CD8 T cells.

# %%
with mpl.rc_context({'figure.figsize': (12, 4)}):
    sc.pl.rank_genes_groups_violin(adata, groups='16', n_genes=20)

# %% [markdown]
# Simple violin plots grouped by cluster also offer a good visualisation:

# %%
sc.pl.violin(adata, ['CST3', 'NKG7', 'PPBP'], groupby='leiden', stripplot = False )

# %% [markdown]
# Dotplot can generate a bird's eye view of expression of all selected marker genes in all clusters.

# %%
sc.pl.dotplot(adata, marker_genes, groupby='leiden')

# %% [markdown]
# Yet another option is given by stacked violin plot.

# %%
sc.pl.stacked_violin(adata, marker_genes, groupby='leiden', rotation=90)
# %% 
adata.write("pbmc10k_scanpy.h5ad")

# %% [markdown]
# -------
# %% [markdown]
# ### PART 5. Cell ontology-based automatic cell type prediction with `CellO`

# %% [markdown]

# %%
#! sudo apt-get -y install graphviz
! sudo apt-get -y install libgraphviz-dev
#! pip install --no-input -U cello-classify --no-cache

# %% [markdown]
# After this (we assume we already imported `pandas`, `numpy`, `scanpy`, and `AnnData` - see above), we run

# %%
import os
import cello

# %% [markdown]
# First of all, let's "recall" the raw data to adata again; we'll need the lists of retained cells and genes for the downstream analysis, but we'll do the normalization and variable gene selection from the scratch. 

# %%
adata = adata.raw.to_adata()
adata

# %%
#new_adata = sc.read_10x_mtx('../data/soupX_pbmc10k_filt', cache = True)
new_adata = adata.copy()

# %%
new_adata = new_adata[adata.obs_names,adata.var_names]
new_adata

# %% [markdown]
# Let's normalize the data and find the highly variable genes - note the unusually high number of HVGs (10k). This follows the suggestion in the `CellO` tutorial:

# %%
sc.pp.normalize_total(new_adata, target_sum=1e6)
sc.pp.log1p(new_adata)
sc.pp.highly_variable_genes(new_adata, n_top_genes=10000)
new_adata

# %% [markdown]
# Now let's cluster the data. Over-clustering is recommended here, for the reasons that would become clear below:

# %%
sc.pp.pca(new_adata, n_comps=50, use_highly_variable=True)
sc.pp.neighbors(new_adata, n_neighbors=15)
sc.tl.leiden(new_adata, resolution=2.0)

# %% [markdown]
# Currently, `CellO` only supports dense matrix format, so we need to do the transformation below. Note that it needs to be done after all `Scanpy` operations.

# %%
new_adata.X = new_adata.X.todense()

# %% [markdown]
# Finally, let's set the `CellO` resource directory to the dir where this notebook is located. NOTE: You will need at least 5 Gb of free disk space for the following code to work. 

# %%
cello_resource_loc = os.getcwd()

# %% [markdown]
# Following this, we can train the model using `CellO` using the following command. This runs for about an hour, and will make a file called `10k_pbmc.model.dill`.
# 
# ```
# model_prefix = "10k_pbmc"
# 
# cello.scanpy_cello(
#     new_adata, clust_key='leiden', rsrc_loc=cello_resource_loc, 
#     out_prefix=model_prefix, log_dir=os.getcwd()
# )
# ```

# %% [markdown]
# Instead, we will download a pre-made model to save time: 

# %%
! wget https://www.dropbox.com/s/x8n6jxa9ygw3s1w/10k_pbmc.model.dill

# %% [markdown]
# Following this, we can run the cell type annotation command with a pre-made model: 

# %%
model_prefix = "10k_pbmc"

cello.scanpy_cello(
    new_adata, 
    clust_key='leiden',
    rsrc_loc=cello_resource_loc,  
    model_file=f'{model_prefix}.model.dill'
)

# %% [markdown]
# Let's run UMAP and visualized the obtained clusters. Quite clearly this is too many, but it's OK: 

# %%
sc.tl.umap(new_adata)

# %%
fig = sc.pl.umap(new_adata, color='leiden', title='Clusters', return_fig=True)

# %% [markdown]
# The following example can be useful to save the plots in the pdf format: 

# %%
out_file = '10k_pbmc_clusters_res2.0.pdf' # <-- Name of the output file

fig.savefig(out_file, bbox_inches='tight', format='pdf')

# %% [markdown]
# Let's visualize a most specific cell type for each cluster: 

# %%
import matplotlib as mpl

with mpl.rc_context({'figure.figsize': (10, 10)}):
    sc.pl.umap(new_adata, color='Most specific cell type',palette="Paired")

# %% [markdown]
# We can see that some of these assignments are wrong; for example, what's marked as "dendiritic cell" in the plot above is actually a population of non-classical monocytes (as shown by expression of FCGR3A). At the same time, what's labelled as "mononuclear cell" is in fact two populations of dendritic cells - myeloid DCs (expressing FCER1A) and plasmacytoid DCs (expressing LILRA4):

# %%
sc.pl.umap(new_adata, color=["FCGR3A","FCER1A","LILRA4"])

# %% [markdown]
# Assigned probability can be used to visualize the probability of cell being a certain cell type:  

# %%
with mpl.rc_context({'figure.figsize': (10, 10)}):
    sc.pl.umap(new_adata, color='natural killer cell (probability)', vmin=0.0, vmax=1.0, return_fig=True)

# %% [markdown]
# We can also use more abstract ontologies, e.g. to identify nucleate cells (platelets don't have nuclei!):

# %%
with mpl.rc_context({'figure.figsize': (10, 10)}):
    sc.pl.umap(new_adata, color='nucleate cell (probability)', vmin=0.0, vmax=1.0, return_fig=True)

# %% [markdown]
# Finally, we can use the result of binary classification:

# %%
with mpl.rc_context({'figure.figsize': (10, 10)}):
    sc.pl.umap(new_adata, color='natural killer cell (binary)', vmin=0.0, vmax=1.0, return_fig=True)

# %% [markdown]
# Additionally, we can visualize cell type probabilities assigned to a specific cluster overlaid on the Cell Ontology graph:

# %%
fig, ax = cello.cello_probs(new_adata, '11', cello_resource_loc, 0.5, clust_key='leiden');




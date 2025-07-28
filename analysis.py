import numpy as np
import pandas as pd
import scanpy as sc

# lines 6 - 29 from: https://www.danli.org/2021/02/03/single-cell-data-analysis-using-scanpy/ 
sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

adata = sc.read_10x_mtx(
    "datasets", 
    var_names="gene_symbols",
    cache=True
)

adata.var_names_make_unique()

#filter 20 highest expressed genes
#sc.pl.highest_expr_genes(adata, n_top=20, )

#filtering out the cells with low gene expression/genes that don't show up in many cells
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

#annotate mitochondrial genes as 'mt' and calculate qc metrics based on 
adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
adata.var['ribo'] = adata.var_names.str.startswith(("RPS", "RPL"))
adata.var['hb'] = adata.var_names.str.startswith("^HB[^(P)]") #regular expression that chooses all leters after HB besides the capital letter P
#regex is used because pseudogenes have names like HBP1, HBP2, etc 

sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], percent_top=None, log1p=False, inplace=True)
# sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
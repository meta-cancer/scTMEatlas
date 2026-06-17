# Spatial scatter visualization of cell-type annotations
import squidpy as sq

sq.pl.spatial_scatter(
    adata,
    library_id='spatial',
    shape=None,
    color=['Cell_type']
)

# Non‑parametric differential expression with violin plots
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

for f in features:

    p = mannwhitneyu(
        adata.obs.loc[adata.obs['Neu16_status']=='Proximal to Neu16', f],
        adata.obs.loc[adata.obs['Neu16_status']=='Distal to Neu16', f]
    ).pvalue

    sns.violinplot(
        data=adata.obs,
        x='Neu16_status',
        y=f,
        inner='box'
    )

    plt.title(f'{'***' if p<0.001 else '**' if p<0.01 else '*' if p<0.05 else 'ns'}')
    plt.show()

# Tangram cell‑cluster projection onto spatial coordinates
import tangram as tg

tg.plot_cell_annotation_sc(
    adata,
    ['Cell_clusters'],
    perc=0.03,
)


# Clustered heatmap of the spatial cellular modules
g = sns.clustermap(
    spatial_matrix,
    cmap='RdBu_r',
    row_colors=row_colors,
    col_colors=col_colors,
    row_cluster=False,
    col_cluster=False,
    figsize=(20, 20)
)



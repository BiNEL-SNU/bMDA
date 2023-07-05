# %%
import os
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
from PIL import Image



# %% read dataframe
order_df = pd.read_csv("data/slide31 phylogeny.csv", dtype = str).dropna(subset={"markingNo"}).rename(columns = {"clusterPhyloSNV_WESorder" : "order"})[["name", "markingNo", "order"]]

position_df = pd.read_csv("data/31 - 190422_H&E_phylo_custom.txt", sep = "\t").\
    rename(columns = {"X1" : "X1_HNE",
                      "Y1" : "Y1_HNE",
                      "X2" : "X2_HNE",
                      "Y2" : "Y2_HNE",
                      "NOTE" : "NOTE_HNE"})
position_df["markingNo"] = position_df.apply(lambda x : x["NOTE_HNE"].split("_", 1)[0], axis = 1)

spatial_df = position_df.merge(order_df, on = "markingNo").astype({"order": float}).sort_values(by = "order").dropna(subset = ["order"])

for i, r in spatial_df.iterrows():
    p1 = (r["X1_HNE"], r["Y1_HNE"])
    p2 = (r["X2_HNE"], r["Y2_HNE"])
    c = ((p1[0] + p2[0]) /2, (p1[1] + p2[1]) /2)
    spatial_df.loc[i, 'x'] = int(c[0])
    spatial_df.loc[i, 'y'] = int(c[1])
print(spatial_df)

snv_data = pd.read_csv("data/190422 WES SNVs.csv")[spatial_df['name']].transpose()


# %% initialize adata
print(snv_data)
snv_np = snv_data.fillna(0).to_numpy()
spatial_np = spatial_df.set_index("name", drop = True)[["x", "y"]].to_numpy()

## move by offset
spatial_np = spatial_np - [1680, 1180]


# %% add artifitial root
snv_np =  np.concatenate([np.zeros((1, snv_np.shape[1])), snv_np], axis=0)
spatial_np =  np.concatenate([np.zeros((1, spatial_np.shape[1])), spatial_np], axis=0)

adata = sc.AnnData(snv_np, obsm={"spatial": spatial_np})
adata.uns['iroot'] = 0

# %%
sc.pp.pca(adata)

sc.pp.neighbors(adata)
sc.tl.leiden(adata)

sc.tl.paga(adata)

sc.tl.diffmap(adata)
sc.tl.dpt(adata)

# %% graph plots
sc.pl.paga(adata, color=['leiden'])

sc.tl.draw_graph(adata, init_pos='paga')
sc.pl.draw_graph(adata, color=['leiden'], legend_loc='on data')

sc.pl.paga_compare(
    adata, threshold=0.03, title='', right_margin=0.2, size=10, edge_width_scale=0.5,
    legend_fontsize=12, fontsize=12, frameon=False, edges=True, save=True)
# %%
pc_c = sc.metrics.morans_i(adata, obsm="X_pca")


# %% plot by weighted average of directions
root = np.argmin(adata.obs["dpt_pseudotime"])

# Get the order of cells based on pseudo-time
order = np.argsort(adata.obs["dpt_pseudotime"])


weighted_directions = []

# Calculate weighted average direction vectors
window_size = 3  # Adjust the window size as needed
# for i in range(window_size, len(order)):
for i in range(len(order)-1):
    start = i
    end = min(i + 1 + window_size, len(order))
    directions = spatial_np[order[start+1:end]] - spatial_np[order[start]]
    weights = np.abs(adata.obs["dpt_pseudotime"].values[order[start+1:end]] - adata.obs["dpt_pseudotime"].values[order[start]])
    weighted_direction = np.average(directions, axis=0, weights=weights)
    weighted_directions.append(weighted_direction)

# Normalize weighted average direction vectors
weighted_directions = np.array(weighted_directions)
# weighted_directions /= np.linalg.norm(weighted_directions, axis=1, keepdims=True)
weighted_directions = weighted_directions * 0.13


color_map = ["black"] + spatial_df['color'].tolist()


# Plot the nodes and directions
fig, ax = plt.subplots(figsize=(8, 8), frameon=False)
image = Image.open(r"data/31 - 190422_contour_targetRegion.jpg")#.convert("RGB")
ax.imshow(image.convert('L'), cmap='gray', vmin=0, vmax=255)
ax.axis('off')



# Plot nodes (spatial locations)
radius = 100
ax.scatter(spatial_np[1:, 0], spatial_np[1:, 1], color=color_map[1:], s=radius)

# Plot directions as arrows
for i, direction in enumerate(weighted_directions):
    if i == 0:
        continue
    start = spatial_np[order[i]]
    end = start + direction

    center_x, center_y = start
    dx, dy = direction[0], direction[1]
    arrow_angle_rad = np.arctan2(dy, dx)
    start_x = center_x + radius * np.cos(arrow_angle_rad) * 0.75 #move arrow to start from outside of circle
    start_y = center_y + radius * np.sin(arrow_angle_rad) * 0.75
    ax.arrow(start_x, start_y, dx, dy, color='black', width=20, head_width = 100, head_length=75)


# Set plot title and axis labels
# ax.set_title('Evolution Graph')
# ax.set_xlabel('X')
# ax.set_ylabel('Y')

if not os.path.exists("results"):
    os.makedirs("results")

plt.savefig(r'results/output.pdf', format='pdf')

plt.show()



# %%

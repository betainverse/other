#! /usr/bin/env python
from pandas import Series, DataFrame
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram

# Read in all the data

print "reading CstR"
df_cstr = pd.read_csv('dCstR_NWMN.csv',usecols=['locus','log2FoldChange','padj'])
print "reading HS"
df_hs = pd.read_csv('HS_NWMN.csv',usecols=['locus','log2FoldChange','padj'])
print "reading CP"
df_cp = pd.read_csv('CP_NWMN.csv',usecols=['locus','log2FoldChange','padj'])



#solve the colorbar problem
import matplotlib
from mpl_toolkits.axes_grid1 import AxesGrid

def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormaps dynamic range to be at zero

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormaps range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highets point in the colormaps range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap

shifted_bwr = shiftedColorMap(matplotlib.cm.bwr,midpoint=0.4,stop=1,name='shifted')



#merge tables according to locus tag, keeping only padj < 0.05 and log2FoldChange >=2.

newcols = []
for col in df_cp.columns:
    if 'locus' not in col:
        col= col+'_cp'
    newcols.append(col)
df_cp.columns = newcols
df_cstr_hs = pd.merge(df_cstr,df_hs, on='locus',suffixes=('_cstr','_hs'))
#print df_cstr_hs[:2]
df_all = pd.merge(df_cstr_hs,df_cp, on='locus')
#print df_all[:2]

significant = (df_all['padj_cstr']<0.05) & (df_all['padj_hs']<0.05) & (df_all['padj_cp']<0.05)
changed = (np.abs(df_all['log2FoldChange_cstr'])>1) & (np.abs(df_all['log2FoldChange_hs'])>1) & (np.abs(df_all['log2FoldChange_cp'])>1)
#changed = (np.abs(df_all['log2FoldChange_cstr'])>1)
#significant = (df_all['padj_cstr']<0.05)
df_padj = df_all[significant & changed]
#print df_padj
#print df_padj.shape
df_change = df_padj[['locus','log2FoldChange_cstr','log2FoldChange_hs','log2FoldChange_cp']]
df_change['locus'] = df_change['locus'].apply(lambda x: int(str(x)[-4:]))
df_change = df_change.set_index('locus')
print df_change
print df_change.shape
loci = list(df_change.index)
print loci
from scipy.spatial.distance import pdist,squareform

row_dist = pd.DataFrame(squareform(pdist(df_change, metric='euclidean')), columns=loci, index=loci)
print row_dist

from scipy.cluster.hierarchy import linkage

row_clusters = linkage(pdist(df_change, metric='euclidean'), method='complete')
#pd.DataFrame(row_clusters, 
#             columns=['row label 1', 'row label 2', 'distance', 'no. of items in clust.'],
#             index=['cluster %d' %(i+1) for i in range(row_clusters.shape[0])])
print row_clusters


row_dendr = dendrogram(row_clusters, labels=list(loci))

fig = plt.figure()

ax = fig.add_subplot(111)

cax = ax.matshow(df_change, interpolation='nearest',cmap='bwr')
fig.colorbar(cax)

print list(loci)
ax.set_xticklabels([''] + ['dCstR','HS','CP'])
ax.set_yticks(range(27))
ax.set_yticklabels(list(loci))

plt.savefig('heatmap.pdf')
plt.show()


from scipy.cluster import hierarchy
# makes dendrogram black (1)
hierarchy.set_link_color_palette(['black'])

# plot row dendrogram
fig = plt.figure(figsize=(8,8))
axd = fig.add_axes([0.09,0.1,0.2,0.6])
row_dendr = dendrogram(row_clusters, orientation='right') # makes dendrogram black (2))

# reorder data with respect to clustering
#df_rowclust = df_change[['log2FoldChange_cstr','log2FoldChange_hs','log2FoldChange_cp']].ix[row_dendr['leaves'][::-1]]
ordered_loci = [list(loci)[x] for x in row_dendr['leaves']]
df_rowclust = df_change.ix[ordered_loci]
print row_dendr['leaves']

print ordered_loci
print df_rowclust
axd.set_xticks([])
#axd.set_yticks([])

# remove axes spines from dendrogram
for i in axd.spines.values():
        i.set_visible(False)

# reorder rows with respect to the clustering
#df_rowclust = df_change.ix[ordered_loci]
#print df_rowclust
# plot heatmap
axm = fig.add_axes([0.26,0.1,0.6,0.6]) # x-pos, y-pos, width, height
cax = axm.matshow(df_rowclust, cmap=shifted_bwr)
fig.colorbar(cax)
#axm.set_xticklabels([''] + list(df_rowclust.columns))
axm.set_xticklabels([''] + ['dCstR','HS','CP'])
axm.set_yticks(range(27))
axm.set_yticklabels(ordered_loci)

plt.savefig('dendr.pdf')
plt.show()


#!/usr/bin/env python
"""
Read in a table of processed RNAseq output and a gbk genome sequence file.
Search for clusters of significantly up- or downregulated genes.
Use the annotated genome file to draw the region of the genome containing
the cluster.
Color genes according to significance or degree of up- or downregulation.
"""

from pandas import DataFrame
import pandas as pd
#import numpy as np
from reportlab.lib import colors
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from reportlab.lib.pagesizes import letter
from math import ceil


max_p = 0.01
intermed_p = 0.0000001 # Close to the least significant of the cst operon (HS, not including tauE or CstR
#csvpath = '/Users/edmonds/git/other/arrows/HS_vs_Control_proteins.csv'
#csvpath = '/Users/edmonds/git/other/arrows/4CsR_vs_Control_proteins.csv'
#outputbase = '4CstR'
csvpath = '/Users/edmonds/git/other/arrows/Nit_vs_Control_proteins.csv'
outputbase = 'Nit'
genomepath = '/Users/edmonds/git/other/arrows/AP009351.gb'
    
#def get_gene_identifier(df_row):
#    return df_row.name

def parsetable(filepath):
    return pd.read_csv(filepath)[['Gene','log2FoldChange','pvalue','begin','end','strand','locus','description']]

#def getneighbors(index,allvalues):
#    neighbors = [index]
#    return 0

def expandbroad(thesevalues,allvalues):
    first = thesevalues[0]
    last = thesevalues[len(thesevalues)-1]
    if first-1 in allvalues:
        return expandbroad([first-1]+thesevalues,allvalues)
    elif last+1 in allvalues:
        return expandbroad(thesevalues+[last+1],allvalues)
    elif first-2 in allvalues:
        return expandbroad([first-2,first-1]+thesevalues,allvalues)
    elif last+2 in allvalues:
        return expandbroad(thesevalues+[last+1,last+2],allvalues)
    elif first-3 in allvalues:
        return expandbroad([first-3,first-2,first-1]+thesevalues,allvalues)
    elif last+3 in allvalues:
        return expandbroad(thesevalues+[last+1,last+2,last+3],allvalues)
    elif first-4 in allvalues:
        return expandbroad([first-4,first-3,first-2,first-1]+thesevalues,allvalues)
    elif last+4 in allvalues:
        return expandbroad(thesevalues+[last+1,last+2,last+3],allvalues)
    else:
        return thesevalues

def expandnarrow(thesevalues,allvalues):
    first = thesevalues[0]
    last = thesevalues[len(thesevalues)-1]
    if first-1 in allvalues:
        return expandbroad([first-1]+thesevalues,allvalues)
    elif last+1 in allvalues:
        return expandbroad(thesevalues+[last+1],allvalues)
    elif first-2 in allvalues:
        return expandbroad([first-2,first-1]+thesevalues,allvalues)
    elif last+2 in allvalues:
        return expandbroad(thesevalues+[last+1,last+2],allvalues)
    elif first-3 in allvalues:
        return expandbroad([first-3,first-2,first-1]+thesevalues,allvalues)
    elif last+3 in allvalues:
        return expandbroad(thesevalues+[last+1,last+2,last+3],allvalues)
    else:
        return thesevalues

def findclusters(alldata):
    sig = alldata[alldata.pvalue < max_p]
    highsig = sig[sig.pvalue < intermed_p]
    sigvals = sig.index.values
    highsigvals = highsig.index.values
    clusters = []
    for value in sigvals:
        if value in highsigvals:
            cluster = expandbroad([value],sigvals)
        else:
            cluster = expandnarrow([value],sigvals)
        if (len(cluster)>1 or cluster[0] in highsigvals) and cluster not in clusters:
            clusters.append(cluster)
    return clusters

def padclusters(clusters):
    newclusters = []
    for cluster in clusters:
        first = cluster[0]
        last = cluster[len(cluster)-1]
        newcluster = [first-1]+cluster+[last+1]
        newclusters.append(newcluster)
    return newclusters

def cluster2ids(cluster,alldata):
    return [alldata.iloc[x].Gene for x in cluster]



## def printclusters(clusters,alldata):
##     genome = SeqIO.read(genomepath,"genbank")
##     gd_diagram = GenomeDiagram.Diagram(genome.id)
##     sizes = [alldata.iloc[cluster[len(cluster)-1]].end -
##              alldata.iloc[cluster[0]].begin for cluster in clusters]
##     max_size = max(sizes) + 200
##     for cluster in clusters:
##         gd_track_for_features = gd_diagram.new_track(i,name=str(cluster[0]),
##                                                      height=0.5)
##         drawcluster(genome,cluster,alldata,gd_track_for_features)
##     gd_diagram.draw(format="linear", pagesize=letter, fragments = 1,
##                     start=0)
##     gd_diagram.write("out.pdf","PDF")

def printclusters(clusters,alldata):
    genome = SeqIO.read(genomepath,"genbank")
    sizes = [alldata.iloc[cluster[len(cluster)-1]].end -
             alldata.iloc[cluster[0]].begin for cluster in clusters]
    max_size = max(sizes) + 200
    print "max size:", max_size
    page_num = 1
    for group in groupclusters(clusters,15):
        printpage(group,alldata,genome,maxsize = max_size,
                  pagenum=page_num,filenamebase = outputbase)
        page_num = page_num+1

def groupclusters(biglist,maxperpage):
    numpages = ceil(float(len(biglist))/maxperpage)
    pages = []
    newlist = []
    print len(biglist)
    for x in biglist:
        if len(newlist) >= maxperpage:
            pages.append(newlist)
            newlist = []
        newlist.append(x)
    pages.append(newlist)
    return pages

def printpage(clusters,alldata,genome,maxsize=0,pagenum=0,filenamebase = 'out'):
    gd_diagram = GenomeDiagram.Diagram(genome.id)
    if maxsize == 0:
        sizes = [alldata.iloc[cluster[len(cluster)-1]].end -
                 alldata.iloc[cluster[0]].begin for cluster in clusters]
        maxsize = max(sizes) + 200
    i=1
    for cluster in clusters:
        gd_track_for_features = gd_diagram.new_track(i,name=str(cluster[0]))
        drawcluster(genome,cluster,alldata,gd_track_for_features)
        i=i+1
    while i<15:
        gd_track_for_features = gd_diagram.new_track(i,name=str(cluster[0]))
        i=i+1
    gd_diagram.draw(format="linear", pagesize=letter, fragments = 1,
                    start=0,end=maxsize,orientation='portrait')
    gd_diagram.write("%s%d.pdf"%(filenamebase,pagenum),"PDF")

        
def drawcluster(genome,cluster,alldata,track):
    cluster_ids = cluster2ids(cluster,alldata)
    #print cluster_ids
    index = dict(zip(cluster_ids, cluster))
    cluster_start = alldata.iloc[cluster[0]].begin-100
    cluster_end = alldata.iloc[cluster[len(cluster)-1]].end+100
    region = genome[cluster_start:cluster_end]
    gd_feature_set = track.new_set()
    for feature in region.features:
        if 'protein_id' in feature.qualifiers.keys():
            thisid = feature.qualifiers['protein_id'][0]
            #print alldata.iloc[index[thisid]].pvalue
            if thisid in cluster_ids:
                thiscolor=colors.white
                pvalue = alldata.loc[index[thisid]].pvalue
                log2foldchange = alldata.loc[index[thisid]].log2FoldChange
                if pvalue < intermed_p:
                    if log2foldchange < 0:
                        thiscolor = colors.blue
                    else:
                        thiscolor = colors.red
                elif pvalue < max_p:
                    if log2foldchange < 0:
                        thiscolor = colors.lightblue
                    else:
                        thiscolor = colors.pink
                else:
                    thiscolor = colors.white
                gd_feature_set.add_feature(feature,sigil="BIGARROW",
                                           color=thiscolor, label=True,
                                           label_size = 4, label_angle=0,
                                           arrowhead_length=0.25,
                                           arrowshaft_height=0.4)


## def printcluster(cluster,alldata):
##     cluster_ids = cluster2ids(cluster,alldata)
##     print cluster_ids
##     cluster_start = alldata.iloc[cluster[0]].begin-100
##     cluster_end = alldata.iloc[cluster[len(cluster)-1]].end+100
##     #fake_end = (cluster_end - cluster_start)*10
##     genome = SeqIO.read(genomepath,"genbank")
##     region = genome[cluster_start:cluster_end]
##     gd_diagram = GenomeDiagram.Diagram(genome.id)
##     gd_track_for_features = gd_diagram.new_track(1,
##                                                  name="Annotated Features",
##                                                  height=0.5)

##     for i in range(2,6):
##         gd_diagram.new_track(i,name=str(i),height=0.5)
    
##     gd_feature_set = gd_track_for_features.new_set()
##     for feature in region.features:
##         if 'protein_id' in feature.qualifiers.keys():
##             id = feature.qualifiers['protein_id'][0]
##             print id
##             if id in cluster_ids:
##                 gd_feature_set.add_feature(feature,sigil="BIGARROW",
##                                            color=colors.blue, label=True,
##                                            label_size = 6, label_angle=90,
##                                            arrowhead_length=0.5,arrowshaft_height=0.7)
##                 #gd_feature_set.add_feature(feature, sigil="ARROW", color="red",
##                 #           arrowhead_length=1)
##                 print feature

##     gd_diagram.draw(format="linear", pagesize=letter, fragments = 1,
##                     start=0)
##     gd_diagram.write("out.pdf","PDF")

## def printgenes(features):
##     record = SeqIO.read("AP009351.gb", "genbank")

##     gd_diagram = GenomeDiagram.Diagram(record.id)
##     gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
##     gd_feature_set = gd_track_for_features.new_set()
##     genome = SeqIO.read(genomepath,"genbank")
##     for feature in record.features:
##         if feature.type != "gene":
##             #Exclude this feature
##             continue
##         if len(gd_feature_set) % 2 == 0:
##             color = colors.blue
##         else:
##             color = colors.lightblue
##         gd_feature_set.add_feature(feature, sigil="BIGARROW",
##                                color=color, label=True,
##                                label_size = 6, label_angle=90)
##     gd_diagram.draw(format="linear", pagesize=letter, fragments=10,
##                 start=0, end=len(record)/10)
##     gd_diagram.write("90.pdf", "PDF")


def main():
    alldata = parsetable(csvpath)
#    print highsig.iloc[0].name
#    print alldata.iloc[highsig.iloc[0].name]
    clusters = findclusters(alldata)
    paddedclusters = padclusters(clusters)
#    for cluster in paddedclusters:
#        print cluster2ids(cluster,alldata)

    printclusters(paddedclusters,alldata)
    #for group in groupclusters(paddedclusters,10):
    #    for cluster in group:
    #        print cluster
    #    print "new group"
    
main()

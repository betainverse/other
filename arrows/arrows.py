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


max_p = 0.01
intermed_p = 0.0000001 # Close to the least significant of the cst operon (HS, not including tauE or CstR
csvpath = '/Users/edmonds/git/other/arrows/HS_vs_Control_proteins.csv'
genomepath = '/Users/edmonds/git/other/arrows/AP009351.gb'

def printgenes(features):
    genome = SeqIO.read(genomepath,"genbank")
    for feature in record.features:
        if feature.type != "gene":
            #Exclude this feature
            continue
        if len(gd_feature_set) % 2 == 0:
            color = colors.blue
        else:
            color = colors.lightblue
        gd_feature_set.add_feature(feature, sigil="BIGARROW",
                               color=color, label=True,
                               label_size = 6, label_angle=90)
    gd_diagram.draw(format="linear", pagesize=letter, fragments=10,
                start=0, end=len(record)/10)
    gd_diagram.write("90.pdf", "PDF")
    
def get_gene_identifier(df_row):
    return df_row.name

def parsetable(filepath):
    return pd.read_csv(filepath,index_col='Gene')[['log2FoldChange','pvalue','begin','end','strand','locus','description']]



def main():
    alldata = parsetable(csvpath)
    sig = alldata[alldata.pvalue < max_p]
    highsig = sig[sig.pvalue < intermed_p]
    print highsig.ix[0].name
    print alldata.ix[get_gene_identifier(highsig.ix[0])]


main()

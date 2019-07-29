#!/usr/bin/env python2


import sys
import re
import os
import argparse
import numpy as np
import pandas as pd
from scipy.spatial import distance

def Dist(array1, array2, dist):
    if dist == 'braycurtis':
        return distance.braycurtis(array1, array2)
    elif dist == 'correlation':
        return distance.correlation(array1, array2)
    elif dist == 'mahalanobis':
        return distance.mahalanobis(array1, array2)
    elif dist == 'minkowski':
        return distance.minkowski(array1, array2)
    elif dist == 'seuclidean':
        return distance.seuclidean(array1, array2)
    elif dist == 'sqeuclidean':
        return distance.sqeuclidean(array1, array2)


def GetMitrix(file1, file2, dist, control):
    pd1 = pd.read_csv(file1, sep='\t', header=0, index_col=0)
    pd2 = pd.read_csv(file2, sep='\t', header=0, index_col=0)
    if control == 'col':
        pd1 = pd1.T
        pd2 = pd2.T
    pd_out = pd.DataFrame(0, index=pd1.index, columns=pd2.index)
    for index1, row1 in pd1.iterrows():
        for index2, row2 in pd2.iterrows():
            pd_out.loc[index1, index2] = Dist(row1, row2, dist)
    pd_out.to_csv('correlation.txt', sep='\t')
    return pd_out

def Heatmap(dict_pht):
    dict_pheatmap_arg = dict_pht
    for key in dict_pheatmap_arg:
        try:
            m = int(dict_pheatmap_arg[key])
        except:
            if dict_pheatmap_arg[key] in ['FALSE', 'TRUE', 'T', 'F']:
                pass
            else:
                dict_pheatmap_arg[key] = '"{}"'.format(dict_pheatmap_arg[key])
    if 'fontsize' not in dict_pheatmap_arg:
        dict_pheatmap_arg['fontsize'] = '4'
    if 'border' not in dict_pheatmap_arg:
        dict_pheatmap_arg['border'] = 'FALSE'
    if 'scale' not in dict_pheatmap_arg:
        dict_pheatmap_arg['scale'] = '"row"'
    if 'filename' not in dict_pheatmap_arg:
        out = 'heatmap.pdf'
    else:
        out = dict_pheatmap_arg['filename']
    pt_arg = ', '.join([i+'='+dict_pheatmap_arg[i] for i in dict_pheatmap_arg])
    R_code = open('Heatmap.R', 'w')
    R_code.write("library(pheatmap)\n\
                 dataExpr <- read.table('correlation.txt', sep='\\t', row.names=1, header=T, quote="", comment="", check.names=F)\n\
                 pdf('{}')\n\
                 out <- pheatmap(dataExpr, {})\n\
                 dev.off()\n\
                 ".format(out, pt_arg))
    R_code.close()
    os.system('Rscript Heatmap.R')
    os.system('rm Heatmap.R')

def main():
    parser = argparse.ArgumentParser(description="Sample Correlation")
    parser.add_argument('--f1', help='input file1', required=True)
    parser.add_argument('--f2', help='input file2')
    parser.add_argument('--dist', help='input distance method',\
                        choices=['braycurtis', 'correlation', 'mahalanobis',\
                                'minkowski', 'seuclidean', 'sqeuclidean'],\
                                default='correlation')

    parser.add_argument('--orientation', help='input orientation of input',\
                       choices=['row', 'col'], default='row')
    parser.add_argument('--pheatmap', help='input pheatmap parameter')
    argv=vars(parser.parse_args())
    if not argv['f2']:
        f2 = argv['f1']
    else:
        f2 = argv['f2']
    dict_pht = {}
    if argv['pheatmap']:
        list_split = re.split(',', argv['pheatmap'])
        for line in list_split:
            key, value = re.split(':', line)
            dict_pht[key] = value
    pd_out = GetMitrix(argv['f1'], f2, argv['dist'], argv['orientation'])
    Heatmap(dict_pht)


if __name__ =='__main__':
    main()




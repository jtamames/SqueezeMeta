#!/usr/bin/env python

import cv2
import pandas as pd
import numpy as np

def pca_signatures(signature_file):
    contigs = pd.read_csv(signature_file,header=0,index_col=0)
    #add pseudo count
    contigs_idx = contigs.index
    contigs = contigs.as_matrix()
    contigs += 1
    log_contigs = np.log(contigs / contigs.sum(axis=1,keepdims=True))
    print((log_contigs > 0).any())
    df_log_contigs = pd.DataFrame(log_contigs,index=contigs_idx)
    df_log_contigs.to_csv(signature_file+".log")
    pca = cv2.PCACompute(log_contigs,np.mean(log_contigs,axis=0).reshape((-1,1)))
    return pca

if __name__=="__main__":
    import sys
    signature_file = sys.argv[1]
    pca_results = pca_signatures(signature_file)
    print(pca_results)

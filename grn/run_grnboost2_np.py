import os
import numpy as np
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

data_dir = '/home/brad/data2/rstudio/birds/scRNA/devin_combined/finch_cells/grn/export_to_numpy_glut'
expr_fname = os.path.join(data_dir, '1.1_exprMatrix_filtered_t.txt')
tf_fname = os.path.join(data_dir, '1.1_inputTFs.txt')

if __name__ == '__main__':
    # ex_matrix is a numpy ndarray, which has no notion of column names
    ex_matrix = np.genfromtxt(expr_fname, delimiter='\t', skip_header=1)

    # we read the gene names from the first line of the file
    with open(expr_fname) as file:
        gene_names = [gene.strip() for gene in file.readline().split('\t')]

    # sanity check to verify the ndarray's nr of columns equals the length of the gene_names list
    assert ex_matrix.shape[1] == len(gene_names)

    # tf_names is read using a utility function included in Arboreto
    tf_names = load_tf_names(tf_fname)

    network = grnboost2(expression_data=ex_matrix,
                        gene_names=gene_names,  # specify the gene_names
                        tf_names=tf_names)

    network.to_csv('output.tsv', sep='\t', index=False, header=False)
                                                                    

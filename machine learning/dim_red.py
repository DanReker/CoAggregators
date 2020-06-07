#######
# library imports
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
import numpy as np
import pylab as pl
from sklearn.manifold import MDS

# read drugs and convert to Morgan Fingerprints
all_drugs = Chem.SmilesMolSupplier("../data/drugbank5_approved_names_smiles.tsv",delimiter="\t",smilesColumn=1, nameColumn=2, titleLine=True)
fps_all_drugs = [AllChem.GetMorganFingerprintAsBitVect(m,4,nBits=2048) for m in all_drugs if not m is None]

sel_drugs = Chem.SmilesMolSupplier("../data/selected_drugs_smiles.tsv",delimiter="\t",smilesColumn=1, nameColumn=0, titleLine=True)
fps_sel_drugs = [AllChem.GetMorganFingerprintAsBitVect(m,4,nBits=2048) for m in sel_drugs]

# generate tanimoto DISsimilarity matrix
fps = fps_all_drugs + fps_sel_drugs
tanimoto_dist_matrix = np.array([[1.0-DataStructs.TanimotoSimilarity(fp1,fp2) for fp1 in fps] for fp2 in fps])

# generate color scheme, all drugs in blue and selected in red
cols = ["#0056DD" for i in range(len(fps_all_drugs))] + ["#EE011A" for i in range(len(fps_sel_drugs))]

# dimension reduction via MDS using tanimoto dissimilarities
embedding = MDS(n_components=2,dissimilarity='precomputed')
embedding.fit(tanimoto_dist_matrix)

# plot MDS and print to file
pl.scatter(embedding.embedding_[:,0],embedding.embedding_[:,1],c=cols,s=50,alpha=0.7,edgecolors="lightgray")
pl.savefig("mds_drugs.pdf")
pl.close()

#############################################


# read excipients and convert to Morgan Fingerprints
all_exc = Chem.SmilesMolSupplier("../data/gras_iig.tsv",delimiter="\t",smilesColumn=2, nameColumn=0, titleLine=True)
fps_all_exc = [AllChem.GetMorganFingerprintAsBitVect(m,4,nBits=2048) for m in all_exc if not m is None]

sel_exc = Chem.SmilesMolSupplier("../data/selected_excipients_smiles.tsv",delimiter="\t",smilesColumn=1, nameColumn=0, titleLine=True)
fps_sel_exc = [AllChem.GetMorganFingerprintAsBitVect(m,4,nBits=2048) for m in sel_exc]

# generate tanimoto DISsimilarity matrix
fps2 = fps_all_exc + fps_sel_exc
tanimoto_dist_matrix2 = np.array([[1.0-DataStructs.TanimotoSimilarity(fp1,fp2) for fp1 in fps2] for fp2 in fps2])

# generate color scheme, all excipients in blue and selected in red
cols2 = ["#0056DD" for i in range(len(fps_all_exc))] + ["#EE011A" for i in range(len(fps_sel_exc))]

# dimension reduction via MDS using tanimoto dissimilarities
embedding2 = MDS(n_components=2,dissimilarity='precomputed')
embedding2.fit(tanimoto_dist_matrix2)

# plot MDS and print to file
pl.scatter(embedding2.embedding_[:,0],embedding2.embedding_[:,1],c=cols2,s=50,alpha=0.7,edgecolors="lightgray")
pl.savefig("mds_excs.pdf")
pl.close()





#########
# import libraries
import numpy as np
from sklearn.ensemble import RandomForestClassifier as RF
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors
import math


#########
# feature definition

descr = Descriptors._descList
calc = [x[1] for x in descr]

def describe_mol(mol):
	fp = AllChem.GetMorganFingerprintAsBitVect(mol,4,nBits=2048)
	fp_list = []
	fp_list.extend(fp.ToBitString())
	fp_expl = []
	fp_expl = [float(x) for x in fp_list]
	ds_n = []
	for d in calc:
		v = d(mol)
		if v > np.finfo(np.float32).max: 	# postprocess descriptors for freak large values
			ds_n.append(np.finfo(np.float32).max)
		elif math.isnan(v):
			ds_n.append(np.float32(0.0))
		else:
			ds_n.append(np.float32(v))
    
	return fp_expl + list(ds_n)


#########
# read training data

d_smiles = []
d_names = []
infile = open("../data/selected_drugs_smiles.tsv","rb")
_ = infile.next() # skip header
for line in infile:
	d_smiles += [line.split("\t")[1].strip()]
	d_names += [line.split("\t")[0]]


e_smiles = []
e_names = []
infile = open("../data/selected_excipients_smiles.tsv","rb")
_ = infile.next() # skip header
for line in infile:
	e_smiles += [line.split("\t")[1].strip()]
	e_names += [line.split("\t")[0].strip()]


screen_results = np.loadtxt("../data/screening_data.tsv",delimiter='\t',usecols=(2,),skiprows=1)
screen_drugs = np.loadtxt("../data/screening_data.tsv",delimiter='\t',usecols=(0,),skiprows=1,dtype=object)
screen_excs = np.loadtxt("../data/screening_data.tsv",delimiter='\t',usecols=(1,),skiprows=1,dtype=object)


e_x = [describe_mol(Chem.MolFromSmiles(e)) for e in e_smiles]
d_x = [describe_mol(Chem.MolFromSmiles(d)) for d in d_smiles]

d_x_dict = dict(zip(d_names, d_x))
e_x_dict = dict(zip(e_names, e_x))



x = []
y = []
for j in range(len(screen_results)):
	x += [d_x_dict[screen_drugs[j]] + e_x_dict[screen_excs[j]]]
	y += [screen_results[j]]


x = np.array(x)
y = np.array(y)



#########
# train model
model = RF(n_estimators=500,n_jobs=25)
model.fit(x,y)


#########
# read other drugs and excipients for prediction


drugs = []
drugs_names = []
infile = open('../data/drugbank_selfaggs_smiles.tsv','r')
_ = infile.next() # skip header
for line in infile:
	drugs += [Chem.MolFromSmiles(line.split('\t')[1])]
	drugs_names += [line.split('\t')[0]]



excipients = []
excipients_names = []
infile = open('../data/gras_iig.tsv','r')
_ = infile.next() # skip header
for line in infile:
	if not Chem.MolFromSmiles(line.split('\t')[2]) is None:
		excipients += [Chem.MolFromSmiles(line.split('\t')[2])]
		excipients_names += [line.split('\t')[0]]



other_drugs = []
other_drugs_names = []
infile = open('../data/drugbank5_approved_names_smiles.tsv','r')
_ = infile.next() # skip header
for line in infile:
	if not Chem.MolFromSmiles(line.split('\t')[1]) is None and not line.split("\t")[-1].strip() in d_names + e_names:
		other_drugs += [Chem.MolFromSmiles(line.split('\t')[1])]
		other_drugs_names += [line.split('\t')[2].strip()]


drugs_d = [describe_mol(m) for m in drugs]
excipients_d = [describe_mol(m) for m in excipients]
other_drugs_d = [describe_mol(m) for m in other_drugs]

#########
# run predictions
####
# perform this per drug to avoid memory overflow
# full matrix of drug-excipients might exceed hardware capacity

p1 = []
p2 = []
p = np.array([])
for i in range(len(drugs)):
	x_temp = [] 
	for j in range(len(excipients)):
		x_temp += [drugs_d[i] + excipients_d[j]]
		p1 += [drugs_names[i]]
		p2 += [excipients_names[j]]
    
	p = np.append(p,model.predict_proba(x_temp)[:,1])
	print "yo"


for i in range(len(drugs)):
	x_temp = []
	for j in range(len(other_drugs)):
		x_temp += [drugs_d[i] + other_drugs_d[j]]
		p1 += [drugs_names[i]]
		p2 += [other_drugs_names[j]]
    
	p = np.append(p,model.predict_proba(x_temp)[:,1])
	print "yo"




outfile = open("pair_predictions.tsv","w")
outfile.write("DRUG\tEXCIPIENT\tPREDICTION\n")
for i in range(len(p)):
	if p[i] > 0.2:
		outfile.write( p1[i] + '\t' + p2[i].strip() + "\t" + str(p[i]) + "\n")


outfile.flush()
outfile.close()

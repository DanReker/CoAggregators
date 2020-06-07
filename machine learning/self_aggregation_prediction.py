from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from sklearn.model_selection import KFold
import numpy as np
from openpyxl import load_workbook


# functionaly to describe molecules via fingerprints and descriptors
descr = Descriptors._descList
calc = [x[1] for x in descr]
def describe(mols):
	descrs = []
	for mol in mols:
		fp = AllChem.GetMorganFingerprintAsBitVect(mol,3,nBits=2048)		
		fp_list = []
		fp_list.extend(fp.ToBitString())
		fp_expl = [float(x) for x in fp_list]
		ds_n = []
		for d in calc:
			v = d(mol)
			if v > np.finfo(np.float32).max: 	# postprocess descriptors for freak large values 
				ds_n.append(np.finfo(np.float32).max)
			else:
				ds_n.append(np.float32(v))
		
		descrs += [fp_expl + list(ds_n)];
	
	return descrs




### file available via wget http://bkslab.org/static/files/aggregators/aggregator_hts.xls
### converted to XLSX manually for compatability w OpenPyXL
wb = load_workbook(filename = '../data/aggregator_hts.xlsx')


# process file
smiles = [x.value for x in wb.get_sheet_by_name("Master_List")["A"][1:]]
mols = [Chem.MolFromSmiles(s.encode("utf-8")) for s in smiles]
fps = np.array(describe(mols))
annotation = [x.value for x in wb.get_sheet_by_name("Master_List")["H"][1:]]

def classano(x):
	if x == "AGG":
		return 1
	elif x == "NONAGG":
		return 0
	else:
		return -1

annoclass = np.array([classano(x) for x in annotation]) 

y = annoclass[annoclass != -1]
x = fps[annoclass != -1]


# model evaluation
RF = RandomForestClassifier(n_estimators=1000,n_jobs=1,max_features=None)


def tenfoldcross(x,y,model):
	x = np.array(x)
	y = np.array(y)
	kf = KFold(10,shuffle=True)
	preds = []
	vals  = []
	probs = []
	for train, test in kf.split(x):
		model.fit(x[train],y[train])
		preds = np.append(preds,model.predict(x[test]))
		vals  = np.append(vals,y[test])
		probs = np.append(probs,model.predict_proba(x[test])[:,1])
	
	return [preds,vals,probs]


benchmarks = []
for i in range(10):
	pred_vals, true_vals, pred_prob = tenfoldcross(x,y,RF)
	
	precision, recall, thresholds = metrics.precision_recall_curve(true_vals, pred_prob)
	benchmarks += [[metrics.matthews_corrcoef(true_vals,pred_vals),
	metrics.accuracy_score(true_vals,pred_vals),
	metrics.precision_score(true_vals,pred_vals),
	metrics.f1_score(true_vals,pred_vals),
	metrics.roc_auc_score(true_vals,pred_prob),
	metrics.auc( recall, precision)]]


print "MCC\tAccruacy\tPrecision\tF1\tROCAUC\tAUPR"
print "\t".join([str(f) for f in np.mean(np.array(benchmarks),axis=0)])
print "\t".join([str(f) for f in np.std(np.array(benchmarks),axis=0)])



# fit model on complete dataset
RF.fit(x,y)


# predict new molecules from drugbank
drugbank_mol = []
drugbank_name = []
drugbank_file = open("../data/drugbank5_approved_names_smiles.tsv","r")
 
_ = drugbank_file.next() # skip header
for line in drugbank_file:
	if not Chem.MolFromSmiles(line.split("\t")[1]) is None:
		drugbank_mol += [Chem.MolFromSmiles(line.split("\t")[1])]
		drugbank_name += [line.split("\t")[2].strip()]

drugbank_descr = describe(drugbank_mol)
from sklearn.preprocessing import Imputer
imp = Imputer(missing_values='NaN', strategy='mean', axis=0)
imp.fit(drugbank_descr)
drugbank_descr_impu = imp.transform(drugbank_descr)


predictions_drugbank = RF.predict(drugbank_descr_impu)
probabilities_drugbank = RF.predict_proba(drugbank_descr_impu)[:,1]

outfile = open("../data/drugbank_selfaggs_smiles.tsv","w")
outfile.write("NAME\tSMILES\n")
for i in range(len(drugbank_name)):
	if probabilities_drugbank[i] > 0.2:
		outfile.write(drugbank_name[i] +"\t"+ Chem.MolToSmiles(drugbank_mol[i]) + "\n")

outfile.flush()
outfile.close()

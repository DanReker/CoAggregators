#######
# import libraries 
import numpy as np
import scipy
import scipy.stats
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors
import math


#######
# read data
drug_smiles = []
drug_names = []
infile = open("../data/selected_drugs_smiles.tsv","rb")
_ = infile.next() # skip header
for line in infile:
	drug_smiles += [line.split("\t")[1].strip()]
	drug_names += [line.split("\t")[0]]


exc_smiles = []
exc_names = []
infile = open("../data/selected_excipients_smiles.tsv","rb")
_ = infile.next() # skip header
for line in infile:
	exc_smiles += [line.split("\t")[1].strip()]
	exc_names += [line.split("\t")[0].strip()]


screen_results = np.loadtxt("../data/screening_data.tsv",delimiter='\t',usecols=(2,),skiprows=1)
screen_drugs = np.loadtxt("../data/screening_data.tsv",delimiter='\t',usecols=(0,),skiprows=1,dtype=object)
screen_excs = np.loadtxt("../data/screening_data.tsv",delimiter='\t',usecols=(1,),skiprows=1,dtype=object)


########
# define chemical features for molecular descriptions 
descr = Descriptors._descList
calc = [x[1] for x in descr]
d_name = [x[0] for x in descr]


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



exc_x = [describe_mol(Chem.MolFromSmiles(e)) for e in exc_smiles]
drugs_x = [describe_mol(Chem.MolFromSmiles(d)) for d in drug_smiles]

drug_x_dict = dict(zip(drug_names, drugs_x))
exc_x_dict = dict(zip(exc_names, exc_x))



simulation_results = np.loadtxt("../data/simulation_results.tsv",delimiter='\t',usecols=np.arange(2,21))
simulation_drugs = np.loadtxt("../data/simulation_results.tsv",delimiter='\t',usecols=(0,),dtype=object)
simulation_excs = np.loadtxt("../data/simulation_results.tsv",delimiter='\t',usecols=(1,),dtype=object)

md_dict = {}
for i in range(len(simulation_drugs)):
	if simulation_drugs[i] in md_dict:
		md_dict[simulation_drugs[i]][simulation_excs[i]] = simulation_results[i].tolist()
	else:
		md_dict[simulation_drugs[i]] = {}
		md_dict[simulation_drugs[i]][simulation_excs[i]] = simulation_results[i].tolist()



########
# describe data in terms of chemical features

x = []
y = []
for j in range(len(screen_results)):
	x += [drug_x_dict[screen_drugs[j]] + exc_x_dict[screen_excs[j]] + md_dict[screen_drugs[j]][screen_excs[j]]]
	y += [screen_results[j]]


x = np.array(x)
y = np.array(y)



########
# libraries used for machine learning evaluations
from sklearn.model_selection import KFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.svm import LinearSVC
from sklearn.utils import shuffle
from sklearn.metrics import matthews_corrcoef as mcc
from sklearn.metrics import accuracy_score as ac
from sklearn import metrics
from sklearn.metrics import roc_auc_score as auc

models = [ GaussianNB(),KNeighborsClassifier(3), DecisionTreeClassifier(),MLPClassifier(),LinearSVC()]


########
# function to run 10-fold cross validations
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
		if isinstance(model, LinearSVC):
			probs = np.append(probs,model.decision_function(x[test]))
		else:		
			probs = np.append(probs,model.predict_proba(x[test])[:,1])
	
	return [preds,vals,probs]


########
# function to standardize calculations of benchmark statistics
def evaluate(pred_vals,true_vals,pred_prob):
	precision, recall, thresholds = metrics.precision_recall_curve(true_vals, pred_prob)
	return [mcc(true_vals,pred_vals),
	metrics.f1_score(true_vals,pred_vals),
	metrics.precision_score(true_vals,pred_vals),
	ac(true_vals,pred_vals),
	metrics.roc_auc_score(true_vals,pred_prob),
	metrics.auc( recall, precision)]



#############
# evaluate random forest model

model = RandomForestClassifier(n_estimators=500)
benchmark = []
for rep in range(10):
	true_vals = np.array([])
	pred_vals = np.array([])
	pred_prob = np.array([])
	for drug in drug_names:
		mask = screen_drugs == drug
		x_train = x[np.invert(mask),:]
		y_train = y[np.invert(mask)]
		x_test = x[mask,:]
		y_test = y[mask]
		
		_ = model.fit(x_train,y_train)
		
		true_vals = np.append(true_vals,y_test)
		pred_vals = np.append(pred_vals,model.predict(x_test))
		pred_prob = np.append(pred_prob,model.predict_proba(x_test)[:,1])
		
	
	benchmark += [evaluate(pred_vals,true_vals,pred_prob)]


model = RandomForestClassifier(n_estimators=500)
benchmark_10_cross = []
for i in range(10):
	pred_vals, true_vals, pred_prob = tenfoldcross(x,y,model)
	benchmark_10_cross += [evaluate(pred_vals,true_vals,pred_prob)]	




############
# check performance of model w only simulation or only fingerprints

model = RandomForestClassifier(n_estimators=500)

benchmark_10_cross_noMD = []
for i in range(10):
	pred_vals, true_vals, pred_prob = tenfoldcross(x[:,:-19],y,model)
	benchmark_10_cross_noMD += [evaluate(pred_vals,true_vals,pred_prob)]	

print "\t".join([str(round(temp,2)) for temp in np.mean(benchmark_10_cross,axis=0)])
print "\t".join([str(round(temp,2)) for temp in np.mean(benchmark_10_cross_noMD,axis=0)])
print "\t".join(["{:.0e}".format(scipy.stats.ttest_ind(np.array(benchmark_10_cross)[:,i],np.array(benchmark_10_cross_noMD)[:,i]).pvalue) for i in range(6)])


benchmark_noMD = []
for rep in range(10):
	true_vals = np.array([])
	pred_vals = np.array([])
	pred_prob = np.array([])
	for drug in drug_names:
		mask = screen_drugs == drug
		x_train = x[np.invert(mask),:-19]
		y_train = y[np.invert(mask)]
		x_test = x[mask,:-19]
		y_test = y[mask]
		
		_ = model.fit(x_train,y_train)
		
		true_vals = np.append(true_vals,y_test)
		pred_vals = np.append(pred_vals,model.predict(x_test))
		pred_prob = np.append(pred_prob,model.predict_proba(x_test)[:,1])
		
	
	benchmark_noMD += [evaluate(pred_vals,true_vals,pred_prob)]


print "\t".join([str(round(temp,2)) for temp in np.mean(benchmark,axis=0)])
print "\t".join([str(round(temp,2)) for temp in np.mean(benchmark_noMD,axis=0)])
print "\t".join(["{:.0e}".format(scipy.stats.ttest_ind(np.array(benchmark)[:,i],np.array(benchmark_noMD)[:,i]).pvalue) for i in range(6)])



benchmark_10_cross_onlyMD = []
for i in range(10):
	pred_vals, true_vals, pred_prob = tenfoldcross(x[:,-19:],y,model)
	pred_vals = np.array(pred_prob) > 0.2
	benchmark_10_cross_onlyMD += [evaluate(pred_vals,true_vals,pred_prob)]	


print "\t".join([str(round(temp,2)) for temp in np.mean(benchmark_10_cross_onlyMD,axis=0)])
print "\t".join(["{:.0e}".format(scipy.stats.ttest_ind(np.array(benchmark_10_cross)[:,i],np.array(benchmark_10_cross_onlyMD)[:,i]).pvalue) for i in range(6)])



benchmark_onlyMD = []
for rep in range(10):
	true_vals = np.array([])
	pred_vals = np.array([])
	pred_prob = np.array([])
	for drug in drug_names:
		mask = screen_drugs == drug
		x_train = x[np.invert(mask),-19:]
		y_train = y[np.invert(mask)]
		x_test = x[mask,-19:]
		y_test = y[mask]
		
		_ = model.fit(x_train,y_train)
		
		true_vals = np.append(true_vals,y_test)
		pred_vals = np.append(pred_vals,model.predict(x_test))
		pred_prob = np.append(pred_prob,model.predict_proba(x_test)[:,1])
		
	pred_vals = np.array(pred_prob) > 0.2
	benchmark_onlyMD += [evaluate(pred_vals,true_vals,pred_prob)]

print "\t".join([str(round(temp,2)) for temp in np.mean(benchmark_onlyMD,axis=0)])
print "\t".join(["{:.0e}".format(scipy.stats.ttest_ind(np.array(benchmark)[:,i],np.array(benchmark_onlyMD)[:,i]).pvalue) for i in range(6)])



############
# compare to other models
models = [ GaussianNB(),KNeighborsClassifier(3), DecisionTreeClassifier(),MLPClassifier(),LinearSVC()]

benchmarks_models_10cross = []
for model in models:
	bm_model = []
	for i in range(10):
		pred_vals, true_vals, pred_prob = tenfoldcross(x,y,model)
		bm_model += [evaluate(pred_vals,true_vals,pred_prob)]
	
	benchmarks_models_10cross += [bm_model]
	


for bm in benchmarks_models_10cross:
	print "\t".join([str(round(temp,2)) for temp in np.mean(bm,axis=0)])
	print "\t".join(["{:.0e}".format(scipy.stats.ttest_ind(np.array(benchmark_10_cross)[:,i],np.array(bm)[:,i]).pvalue) for i in range(6)])



benchmarks_models = []
for model in models:
	bm_model = []
	for rep in range(10):
		true_vals = np.array([])
		pred_vals = np.array([])
		pred_prob = np.array([])
		for drug in drug_names:
			mask = screen_drugs == drug
			x_train = x[np.invert(mask),:]
			y_train = y[np.invert(mask)]
			x_test = x[mask,:]
			y_test = y[mask]
			
			_ = model.fit(x_train,y_train)
			
			true_vals = np.append(true_vals,y_test)
			pred_vals = np.append(pred_vals,model.predict(x_test))
			if isinstance(model, LinearSVC):
				pred_prob = np.append(pred_prob,model.decision_function(x_test))
			else:		
				pred_prob = np.append(pred_prob,model.predict_proba(x_test)[:,1])
			
		
		
		bm_model += [evaluate(pred_vals,true_vals,pred_prob)]
	
	benchmarks_models += [bm_model]

	
for bm in benchmarks_models:
	print "\t".join([str(round(temp,2)) for temp in np.mean(bm,axis=0)])
	print "\t".join(["{:.0e}".format(scipy.stats.ttest_ind(np.array(benchmark_10_cross)[:,i],np.array(bm)[:,i]).pvalue) for i in range(6)])






###############
#### adversial shuffle control 

model = RandomForestClassifier(n_estimators=500)
benchmarks_y_shuffle = []
for rep in range(10):
	true_vals = np.array([])
	pred_vals = np.array([])
	pred_prob = np.array([])
	for drug in drug_names:
		mask = screen_drugs == drug
		x_train = x[np.invert(mask),:]
		y_train = np.random.permutation(y[np.invert(mask)])
		x_test = x[mask,:]
		y_test = y[mask]
		
		_ = model.fit(x_train,y_train)
		
		true_vals = np.append(true_vals,y_test)
		pred_vals = np.append(pred_vals,model.predict(x_test))
		pred_prob = np.append(pred_prob,model.predict_proba(x_test)[:,1])
		
	
	benchmarks_y_shuffle += [evaluate(pred_vals,true_vals,pred_prob)]


print "\t".join([str(round(temp,2)) for temp in np.mean(benchmarks_y_shuffle,axis=0)])
print "\t".join(["{:.0e}".format(scipy.stats.ttest_ind(np.array(benchmark)[:,i],np.array(benchmarks_y_shuffle)[:,i]).pvalue) for i in range(6)])


benchmarks_x_shuffle = []
for rep in range(10):
	true_vals = np.array([])
	pred_vals = np.array([])
	pred_prob = np.array([])
	for drug in drug_names:
		mask = screen_drugs == drug
		x_train = x[np.invert(mask),:]
		for x_t in x_train:
			np.random.shuffle(x_t)
		
		y_train = y[np.invert(mask)]
		x_test = x[mask,:]
		y_test = y[mask]
		
		_ = model.fit(x_train,y_train)
		
		true_vals = np.append(true_vals,y_test)
		pred_vals = np.append(pred_vals,model.predict(x_test))
		pred_prob = np.append(pred_prob,model.predict_proba(x_test)[:,1])
		
	
	benchmarks_x_shuffle += [evaluate(pred_vals,true_vals,pred_prob)]


print "\t".join([str(round(temp,2)) for temp in np.mean(benchmarks_x_shuffle,axis=0)])
print "\t".join(["{:.0e}".format(scipy.stats.ttest_ind(np.array(benchmark)[:,i],np.array(benchmarks_x_shuffle)[:,i]).pvalue) for i in range(6)])


benchmark_y_shuffle_10cross = []
for rep in range(10):
	ys = np.random.permutation(y)
	pred_vals, true_vals, pred_prob = tenfoldcross(x,ys,model)
	benchmark_y_shuffle_10cross += [evaluate(pred_vals,true_vals,pred_prob)]

print "\t".join([str(round(temp,2)) for temp in np.mean(benchmark_y_shuffle_10cross,axis=0)])
print "\t".join(["{:.0e}".format(scipy.stats.ttest_ind(np.array(benchmark_10_cross)[:,i],np.array(benchmark_y_shuffle_10cross)[:,i]).pvalue) for i in range(6)])



benchmark_x_shuffle_10cross = []
for rep in range(10):
	xs = x.copy()
	for xs_t in xs:
		np.random.shuffle(xs_t)
	
	pred_vals, true_vals, pred_prob = tenfoldcross(xs,y,model)
	benchmark_x_shuffle_10cross += [evaluate(pred_vals,true_vals,pred_prob)]


print "\t".join([str(round(temp,2)) for temp in np.mean(benchmark_x_shuffle_10cross,axis=0)])
print "\t".join(["{:.0e}".format(scipy.stats.ttest_ind(np.array(benchmark_10_cross)[:,i],np.array(benchmark_x_shuffle_10cross)[:,i]).pvalue) for i in range(6)])



###############
#### feature ablation
from joblib import Parallel, delayed

def ablation(strategy):
	model = RandomForestClassifier(n_estimators=50,oob_score=True)
	assert strategy in ["random", "important"]
	
	x_train = np.array(x).copy()
	y_train = np.array(y).copy()
	mccs = []
	for feat in range(len(x[0])):
		
		_ = model.fit(x_train,y_train)
		
		if strategy == "random":
			x_train = np.delete(x_train,np.random.randint(len(x_train[0])),axis=1)
		elif strategy == "important":
			x_train = np.delete(x_train,np.argmax(model.feature_importances_),axis=1)
		else:
			continue
		
		mccs += [mcc(np.argmax(model.oob_decision_function_,axis=1),y)]
		
	return mccs


random_mccs = Parallel(n_jobs=10)(delayed(ablation)("random") for i in range(10))
important_mccs = Parallel(n_jobs=10)(delayed(ablation)("important") for i in range(10))

np.savetxt("random_ablation_mcc.txt",random_mccs,delimiter='\t')
np.savetxt("important_ablation_mcc.txt",important_mccs,delimiter='\t')




###############
#### feature importance

feature_names = []
feature_names += ["drug_fp_" + str(i) for i in range(2048)]
feature_names += ["drug_" + d for d in d_name]
feature_names += ["excipient_fp_" + str(i) for i in range(2048)]
feature_names += ["excipient_" + d for d in d_name]
feature_names += ["simulation_contact_av", "simulation_kineticE", "simulation_potentialE", "simulation_AvDist_1", "simulation_AvDist_2", "simulation_AvDist_3", "simulation_AvDist_4", "simulation_d(distance)_max_1", "simulation_d(distance)_max_2", "simulation_d(distance)_max_3", "simulation_d(distance)_max_4", "simulation_d(distance)_avg_1", "simulation_d(distance)_avg_2", "simulation_d(distance)_avg_3", "simulation_d(distance)_avg_4", "simulation_Variance_1", "simulation_Variance_2", "simulation_Variance_3", "simulation_Variance_4"]

model = RandomForestClassifier(n_estimators=500)
_ = model.fit(x,y)

for i in range(50):
	index = np.argsort(-model.feature_importances_)[i]
	print "\t".join([str(t) for t in [i+1, feature_names[index], round(model.feature_importances_[index],4)]])

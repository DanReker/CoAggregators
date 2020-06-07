# Code to run machine learning experiments.

Dependency
- Python (2.7)
- Python Packages : RDKit, Scikit-learn, openpyxl


#### self_aggregation_prediction.py
evaluates and runs random forest model to predict self-aggregation propensity of approved drug compounds as candidate compounds for co-aggregation platform. This implementation relies on data from the [Shoichet lab](http://www.bkslab.org/) and is a re-implementation of the algorithm developed in `Feng et al. Nat Chem Biol 1(3),146-8, 2005`

#### retrospective_tests.py
large-scale evaluation of co-aggregation model. includes cross validation experiments, feature ablation, adversarial controls, comparison to other machine learning algorithms.

#### predict_pairs.py
predict new co-aggregation pairs from screening data

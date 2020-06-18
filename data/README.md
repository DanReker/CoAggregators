# Data from high-throughput screening and molecular information to run simulations & machine learning 

All data files are tab-separated and contain a header to define the columns.

#### screening_data.tsv	

Screening results from high-throughput co-aggregation screen.
Contains name of drug and excipient tested. Binary code indicates
whether the pair forms nanoparticles (1) or not (0). Last column
corresponds to percentage size reduction of nanoparticles compared
to unformulated drug.

#### selected_drugs_smiles.tsv	

Chemical SMILES structure of the 16 drugs tested in the
high-throughput co-aggregation screen.

#### selected_excipients_smiles.tsv

Chemical SMILES structure of the 90 excipients tested in the
high-throughput co-aggregation screen.

#### drugbank5_approved_names_smiles.tsv	

Approved compounds from [Drugbank 5.0](https://www.drugbank.ca/).
Serves as input for candidate drug selection as well as 
as resource for other excipient structures.

#### drugbank_selfaggs_smiles.tsv	

Candidate drugs that are predicted to be self-aggregating.

#### gras_iig.tsv	

Approved excipients from [FDA.gov](https://www.fda.gov/).
For more information on dataset curation check
 ` Reker et al. Cell Rep 30(11), 3710-6.e4 (2020) `

#### pair_composition.tsv	

Definition of possible pairs to run simulations. 

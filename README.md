# project

Programme for the retrieving of data, curation of data, training of a machine learning algorithm, and predicting 

#Running of the programme: Environment management

It is necessary to create a conda environment with the file requeriments.yml found in the repository using the command:

conda env create -f requeriments.yml

Followed by its activation:

conda activate ML_MHC_Cytokine

In case you need to deactivate the environment:

conda deactivate

And in case you need to eliminate the environment:

conda remove --name ML_MHC_Cytokine --all

#Requirements for running the programme

It is important that you have a connection to the postgress, therefore, you will need to facilitate to the programme a  PostgreSQL username, password, database and port. 

It is necessary to give to the programme the name of the fasta file with a protein sequence or sequences that you want to use to predict. You can also give the Uniprot Id of the sequence, and the programme will download it. 

It is also necessary to give to the programme on what type of cytokine production you want to conduct the prediction, with IL-10 or IFNg. 



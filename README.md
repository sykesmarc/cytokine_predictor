# project
Programme for the retrieving of data, curation of data, training of a machine learning algorithm, and predicting 

#Running of the programme: Environment management
It is necessary to create a conda environment with the file requeriments.yml found in the repository using the command:
conda env create -f model_environment_gpu.yml
Followed by its activation:
conda activate name

In case you need to deactivate the environment:
conda deactivate
And in case you need to eliminate the environment:
conda remove --name name --all

#Requirements for running the programme
It is important that you have a connection to the postgress, therefore, you will need to facilitate to the programme a  PostgreSQL username, password, database and port. 
It is necessary to give to the programme the name of the file with a protein fasta sequence or sequences that you want to predict if there is binding in MHC and production of cytokines. 



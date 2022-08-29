# project

Programme for the prediction of cytokine release based on IEDB database.  

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

It is important that you have a connection to PostgreSQL, therefore, you will need to facilitate to the programme a  PostgreSQL username, password, database and port. You have to create in advance a database in your PostgreSQL user.  

It is necessary to give to the programme the name of the protein fasta file or sequence ID that you want to use to predict. You can also give the Uniprot Id of the sequence, and the programme will download it. 

It is also necessary to give to the programme on what type of cytokine production you want to conduct the prediction, with IL-10 or IFNg. 



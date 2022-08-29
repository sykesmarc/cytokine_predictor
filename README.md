# Cytokine predictor

Programme for the prediction of cytokine release based on IEDB database.  

#Running of the programme: Environment management

It is necessary to create a conda environment with the file requeriments.yml found in the repository using the command:

conda env create -f requeriments.yml

Followed by its activation:

conda activate ML_Cytokine

In case you need to deactivate the environment:

conda deactivate

And in case you need to eliminate the environment:

conda remove --name ML_Cytokine --all

#Requirements for running the programme

All the information required for the programme is given in the information.txt file. 

It is obligatory to give to the programme a PostgreSQL host, database, user, password and port. You have to create in advance a database in your PostgreSQL user.

It is also obligatory to give an input to the programme, it can be a fasta file or a uniprot ID. 

To specify the cytokine, you need to fill the prediction line with one or more than one cytokine. In order to develop representative models for cytokine prediction. Cytokines with more than 50 samples are: TNFa, IL-2, IL-4, GM-CSF, IL-17A, IL-17, IL-6, IL-5, IL-13 and IFNg. Last revision done September 2022. 

For control usage, you can add the microorganism name in the microorganism line. The control depends on the amount of data present in IEDB. In the case you don't want a control just write no. 

If you want to first predict MHC class II binding of your sequence for different human alleles, selecting the most immunogenic regions of you sample, you have to write yes in the alleles line. 

If you want to predict MHC class II, you need to give a threshold for the cut off in the prediction. As default, you can write 1 and wait for your results. If there are a lot of sequences you can try at 1.3. 

If you don't want to predict any MHC class II binding, write no in the  alleles and threshold lines. 

You can choose the length of the epitopes, with 18 being the biggest length. As default is 15. You have to write the number in the length lines. 



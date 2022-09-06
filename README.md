# Cytokine predictor program

Code repository for the Cytokine release predictor programme. We present a cytokine release predictor programme based solely on epitope's sequence. 

## Abstract 

The program is able to predict cytokine release for the cytokines GM-CSF, IFN-γ, IL-2, IL-4, IL-5, IL-6, IL-10, IL_13, IL-17, IL-17A, TNF-α using already trained models. These models can be updated automatically using the IEDB as the source database. The program can also use a control dataset. In that case, any peptides belonging to the control will be removed from the original training set and the model will be trained again. The program can also make a preliminary prediction of MHC class II binding for the input sequences using the CNN-PepPred tool, in order to identify the potentially immunogenic peptides in the input (if they have not been previously analysed with another tool), before cytokine release prediction. 

## Getting the code

You can get the code using this command:

```
git clone https://github.com/sykesmarc/cytokine_predictor
```

## Environment management

In order to run the program, an environment needs to be set with all the modules required. All modules are found in the requeriments.yml file. 

It is necessary to create a conda environment with the file **_requeriments.yml_** found in the repository using the command:
```
conda env create -f requeriments.yml
```
Followed by its activation:
```
conda activate ML_Cytokine
```
In case you need to deactivate the environment:
```
conda deactivate
```
And in case you need to eliminate the environment:
```
conda remove --name ML_Cytokine --all
```
## Instructions

The user also needs to modify the file __information.txt__ in order to give the desired instructions to the program.

- `update:` Indicates if the models have to be updated (yes) or if the program can use the already built models (no). 
- `host, database, user, password, port:` In the case of a model update, the user must give a PostgreSQL host, database, user, password and port. 
- `input:` Refers to a protein sequence given in a fasta file or uniport ID, which will be used as input for the program. 
- `prediction:` Indicates the program which type of cytokine release prediction assay the user wants to conduct. The cytokine must be written this way: (GM-CSF, IFNg, IL-2, IL-4, IL-5, IL-6, IL-10, IL-13, IL-17A, IL-17, TNFa)
- `microorganism:` Indicates the program if the user wants to add a control dataset to the prediction, if so, the user must give a microorganism name and the update line have to be yes. If the user does not want to add any control has to specify it with no.  
- `alleles:` Indicates if the user wants to make a previous prediction of MHC class II binding before cytokine prediction(yes) or no. 
- `threshold:` When having a previous prediction of MHC class II binding, a threshold for the selection of immunogenic peptides must be given by the user (the most common is 1). If there is no previous prediction of MHC class II binding indicate it with a no. 
- `length:` Indicates the length of epitopes in the prediction of cytokine release. The maximum is 18, and the recommended length is 15.


## Running the programme

In order to run the programme you need to run this command once you have set the environment and modified the _information.txt__ file. 

```
python3 cytokine_predictor.py
```

## Output for cytokine prediction

The output of the program consists of a file named __predictedoutcome.txt__, a csv file with the results of the prediction. This file contains a list of peptides, which come from the input sequences, that have been predicted to produce that particular cytokine, ordered by the predicted outcome. The peptides more likely to induce cytokine production have a higher predicted outcome and appear at the top of the file. This output also gives essential information on the peptides: the protein or sequence they come from, the start and end of the peptide in that protein/sequence, and its characteristic core (reminiscent of the HLA class II binding core). 

## Updating the models

The user can also choose to update the models. To this end, the program downloads the IEDB database, parsers all its xml files and filters the data. After that, it creates a local PostgreSQL table and stores the curated data in it if the models are being generated for the first time, or adds any new data found in the IEDB to the existing PostgreSQL table if the user had already run this process before. Before updating the models for the first time, the user should create a database for its PostgreSQL user. The information required for the PostgreSQL connection is given in the information.txt file.Each updated model will be stored in the folder where the program is running. The performance of the new models can be seen in the cross-validation results provided with the output. 

## Previous prediction of MHC class II binding

The program can also run the CNN-PepPred tool, before cytokine-release prediction, to make a prediction of MHC II binding for the input peptides and for different HLA class II alleles if such analysis has not been done previously with another tool. Then, the program filters the output of each prediction with a given threshold to obtain the peptides predicted to be most immunogenic. These filtered results are subsequently used as input for the prediction of cytokine release.  



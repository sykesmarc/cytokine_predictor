# Cytokine predictor

Code repository for the Cytokine release predictor programme. 

We present a cytokine release predictor programme based solely on epitope's sequence. 

# Abstract 

Programme for the prediction of cytokine release. This programme is able to apply models previsouly created for the prediction cytokine release based solely on epitope's sequence. This programme can also update the models if wanted, apply a previously MHC class II binding prediction, and apply a control dataset. 

# Running the programme: Environment management

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
# Running the programme: information.txt

All the information required for the programme is given in the **_information.txt_** file. 

It is obligatory to give to the programme a PostgreSQL host, database, user, password and port. You have to create in advance a database in your PostgreSQL user.

It is also obligatory to give an input to the programme, it can be a fasta file or a uniprot ID. 

To specify the cytokine, you need to fill the prediction line with one or more than one cytokine. In order to develop representative models for cytokine prediction. Cytokines with more than 50 samples are: TNFa, IL-2, IL-4, GM-CSF, IL-17A, IL-17, IL-6, IL-5, IL-13 and IFNg. Last revision done September 2022. 

For control usage, you can add the microorganism name in the microorganism line. The control depends on the amount of data present in IEDB. In the case you don't want a control just write no. 

If you want to first predict MHC class II binding of your sequence for different human alleles, selecting the most immunogenic regions of you sample, you have to write yes in the alleles line. 

If you want to predict MHC class II, you need to give a threshold for the cut off in the prediction. As default, you can write 1 and wait for your results. If there are a lot of sequences you can try at 1.3. 

If you don't want to predict any MHC class II binding, write no in the  alleles and threshold lines. 

You can choose the length of the epitopes, with 18 being the biggest length. As default is 15. You have to write the number in the length lines. 

# Getting the code

```
git clone https://github.com/sykesmarc/cytokine_predictor
```

# Run the python programme: 

```
python3 cytokine_predictor.py
```


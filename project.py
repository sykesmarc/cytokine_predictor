#Imports
import os
import re
import xml.etree.cElementTree as et
import glob
import shutil
import os
import psycopg2
import sys
from psycopg2.extensions import ISOLATION_LEVEL_AUTOCOMMIT
import subprocess as sp

#Introduction to the programme
print("In order to begin, you need to facilitate the programme the Input and the PostgreSQL connection information")
print("If you have the input in a fasta file write 'FILE'")
print("If you have the Uniprot ID of the input write 'UNIPROT'")
answer = input()
if answer.upper() == "FILE":
    print("Facilitate the name of the file, which should be located where the programme is")
    input1 = input()
else:
    print("Facilitate the Uniprot ID")
    uniprot_ID = input()
    os.system("wget http://www.uniprot.org/uniprot/"+uniprot_ID+".fasta")
    input1 = uniprot_ID+".fasta"

#Obtain information
print("In order to connect to Postgress, wee need you to write this information:")
print("Host")
host1 = input()
print("Database")
database1 = input()
print("User")
user1 = input()
print("Password")
password1 = input()
print("Port")
port1 = input()
print("Thank you so much ;)")

#Prove that it can actually connect
con = psycopg2.connect(
                            host = host1,
                            database = database1,
                            user = user1 ,
                            password = password1,
                            port = port1
                            )
con.set_isolation_level(ISOLATION_LEVEL_AUTOCOMMIT)
con.close()

#IL-10 or IFNg release
print("----------------------------------------------------------------------------")
print("If you want to predict IL-10 production write 'IL-10' and if you want to predict IFNg production write 'IFNg'")
cytokine = input()
print("----------------------------------------------------------------------------")

#Importing importing_files
def importing_files():
    os.system("wget http://immuneepitope.org/doc/iedb_export.zip -O iedb_export.zip")
    os.system("mkdir Database")
    os.system("unzip iedb_export.zip -d Database")
    os.system("wget http://immuneepitope.org/doc/AssayTypeList.zip -O AssayTypeList.zip")
    os.system("unzip AssayTypeList.zip")
    os.system("wget http://immuneepitope.org/doc/OrganismList.zip -O OrganismList.zip")
    os.system("unzip OrganismList.zip")
    os.system("wget http://immuneepitope.org/doc/MhcAlleleNameList.zip -O MhcAlleleNameList.zip")
    os.system("unzip MhcAlleleNameList.zip")
    os.system("rm *.zip")

#Obtain the Dictionaries
def obtain_dicc(file, element1, element2 ):
    dicc = {}
    tree = et.parse(file) #Use the XML parser
    root = tree.getroot() #Build the root of the tree
    for elem in root.iter():
        if elem.tag == element1: #Select the organismId/assayId
            key = elem.text
        if elem.tag == element2: #Select the assay/organism
            value = elem.text
            dicc[key] = value #Add it to a dicctionary
    dicc["NA"] = "NA" #Add a blank option
    return dicc

#Read files
#Read files
def extract_inf(file, assay_dicc, org_dicc, g ,mhc_dicc,f1):
    f = open(file,"rt",encoding="utf-8")
    while True:
        text = f.readline()
        #Matches
        epit = re.search(r'<Epitope>',text)
        epitname = re.search(r'<EpitopeName>',text)
        linseq = re.search(r'<LinearSequence>',text)
        spos = re.search(r'<StartingPosition>',text)
        epos = re.search(r'<EndingPosition>',text)
        genbid = re.search(r'<GenBankId>',text)
        sourceorganismid = re.search(r'<SourceOrganismId>',text)
        assays = re.search(r'<Assays>',text)
        tcellid = re.search(r'<TCellId>',text)
        bcell = re.search(r'<BCell>', text)
        mhcbinding = re.search(r'<MhcBinding>',text)
        organismid = re.search(r'<OrganismId>',text)
        effectorcells = re.search(r'<EffectorCells>',text)
        antigenpresentingcells = re.search(r'<AntigenPresentingCells>',text)
        celltype = re.search(r'<CellType>',text)
        mhcalleleid = re.search(r'<MhcAlleleId>',text)
        assaytypeid = re.search(r'<AssayTypeId>',text)
        qualmeasurements = re.search(r'<QualitativeMeasurement>',text)
        quanmeasurements = re.search(r'<QuantitativeMeasurement>',text)
        finish = re.search(r'</AssayInformation>',text)
        #Read files
        if not text: break
        if epit:#Make sure we are rewritting every epitope we are doing
            epitope = "NA"
            linearsequence = "NA"
            sposition="NA"
            eposition="NA"
            genbankid="NA"
            sourceorganism="NA"
            cellid= "NA"
            hostorganism ="NA"
            ct = "NA"
            ec_celltype = "NA"
            apc_celltype = "NA"
            mhc = "NA"
            mhc_class = "NA"
            assaytype="NA"
            qualitativemeasurements = "NA"
            quantitativemeasurements = "NA"
        #EPITOPE INFORMATION
        if epitname:
            epitope = text.split("EpitopeName>")[1][0:-2]
            epitope = epitope.replace("&lt;","<").replace("&gt;",">").replace(","," ").replace("'"," ").strip("'").strip('"')
        if linseq: linearsequence = str(text.split(">")[1][0:-16])
        if spos: sposition = str(text.split(">")[1][0:-18])
        if epos: eposition = str(text.split(">")[1][0:-16])
        if genbid: genbankid = str(text.split(">")[1][0:-11])
        if sourceorganismid: sourceorganism = str(text.split(">")[1][0:-18])
        #ASSAY INFORMATION
        if assays: #Make sure we are rewriting every epitope
            cellid= "NA"
            hostorganism ="NA"
            ct = "NA"
            ec_celltype = "NA"
            apc_celltype = "NA"
            mhc = "NA"
            mhc_class = "NA"
            assaytype="NA"
            qualitativemeasurements = "NA"
            quantitativemeasurements = "NA"
        if tcellid: cellid = str(text.split(">")[1][0:-9])
        if bcell: cellid = "Bcell"
        if mhcbinding: cellid = "Bcell"
        if organismid: hostorganism = str(text.split(">")[1][0:-12])
        if effectorcells: ct = "EffectorCells"
        if antigenpresentingcells: ct = "AntigenPresentingCells"
        if celltype:
            if ct == "EffectorCells":
                ec_celltype = str(text.split(">")[1][0:-10])
            if ct == "AntigenPresentingCells":
                apc_celltype = str(text.split(">")[1][0:-10])
        if mhcalleleid:
            mhc = str(text.split(">")[1][0:-13])
            if 'HLA-A' in mhc_dicc[mhc] or 'HLA-B' in mhc_dicc[mhc] or 'HLA-C' in mhc_dicc[mhc] or 'HLA-E' in mhc_dicc[mhc] or 'HLA-F' in mhc_dicc[mhc] or 'HLA-G' in mhc_dicc[mhc] or mhc_dicc[mhc] == "HLA class I":
                mhc_class = "MHC I"
            if 'HLA-D' in mhc_dicc[mhc] or mhc_dicc[mhc] == "HLA class II":
                mhc_class = "MHC II"
            if "MHC I" not in mhc_class:
                if "NA" == mhc_dicc[mhc]:
                    mhc_class = "NA"
                else:
                    mhc_class = "Uncaracterized"
        if assaytypeid: assaytype = str(text.split(">")[1][0:-13])
        if qualmeasurements: qualitativemeasurements = str(text.split(">")[1][0:-24])
        if quanmeasurements: quantitativemeasurements = str(text.split(">")[1][0:-25])
        if finish:
            if cellid != "Bcell": #If the experiment is not a B cell experiment, append the row
                if hostorganism == "9606": #Only if it is in humans
                    try:
                        line = "insert into Epitopes (epitopename, linearsequence, startingposition, endingposition, genebankid, sourceorganism, tcellid, hostorganism, eccelltype, apccelltype, mhc, mhctype, assaytype, qualitativemeasurements, quantitativemeasurements) values ("+"'"+ epitope+"' ,"+"'"+linearsequence + "' ,"+ "'"+sposition + "' ," + "'"+eposition+ "' ," + "'"+ genbankid + "' ," +"'"+ org_dicc[sourceorganism] + "' ," +"'"+cellid + "' ," +"'"+ org_dicc[hostorganism]  + "' ," + "'"+ ec_celltype+ "' ," + "'"+ apc_celltype+ "' ," + "'"+ mhc_dicc[mhc]+ "' ," + "'"+ mhc_class + "' ,"+ "'"+ assay_dicc[assaytype] +"' ,"  + "'"+ qualitativemeasurements+ "' ," + "'"+ quantitativemeasurements+"');"
                        f1.write(line+  "\n")
                    except:
                        line = "insert into Epitopes (epitopename, linearsequence, startingposition, endingposition, genebankid, sourceorganism, tcellid, hostorganism, eccelltype, apccelltype, mhc, mhctype, assaytype, qualitativemeasurements, quantitativemeasurements) values ("+"'"+ epitope+"' ,"+"'"+linearsequence + "' ,"+ "'"+sposition + "' ," + "'"+eposition+ "' ," + "'"+ genbankid + "' ," +"'"+ sourceorganism + "' ," +"'"+cellid + "' ," +"'"+ org_dicc[hostorganism] + "' ," + "'"+ ec_celltype+ "' ," + "'"+ apc_celltype+ "' ," + "'"+ mhc_dicc[mhc]+ "' ,"  + "'"+ mhc_class + "' ," + "'"+ assay_dicc[assaytype]+"' ,"  + "'"+ qualitativemeasurements+ "' ," + "'"+ quantitativemeasurements+"');"
                        f1.write(line+  "\n")
    f.close()

#Filter files with T cell experiments
def obtain_file_information(file, assay_dicc, org_dicc,g,mhc_dicc,f1):
    y = 0
    z = 0
    f = open(file, 'rt',encoding="utf-8")
    for i in f.readlines(): #Check if there is any T cell
        match1 = re.search(r'<TCell>',i)
        match2 = re.search(r'<OrganismId>9606</OrganismId>',i)
        if match1: y += 1
        if match2: z += 1
        if y == 1 and z == 1: break
    f.close()
    if y == 0 or z == 0: #If not T cell or human organism, remove file
        os.system("rm -f "+file)
    else: #If T cell, cool
        print("It's a match")
        extract_inf(file, assay_dicc, org_dicc, g ,mhc_dicc,f1)

#Iterate over all the files
def iterate(assay_dicc, org_dicc,mhc_dicc):
    f1 = open("dumping.txt", 'wt',encoding="utf-8" )
    g = 0 #Count to name all the files
    for file in list(glob.glob('Database/*.xml')): #Iterate over all the files
        print(file.strip("Database/"))
        g += 1
        obtain_file_information(file, assay_dicc, org_dicc,g,mhc_dicc,f1)
    f1.close()

#Remove repeated lines
def remove_lines():
    lines = open('dumping.txt', 'r', encoding="utf-8").readlines()
    lines_set = set(lines)
    out = open('database.sql','w', encoding="utf-8")
    for line in lines_set:
            out.write(line)

#Connect to Postgress
def enter_psql(host1, database1, user1, password1, port1): #Connect to postgress
    con = psycopg2.connect(
                            host = host1,
                            database = database1,
                            user = user1 ,
                            password = password1,
                            port = port1
                            )
    con.set_isolation_level(ISOLATION_LEVEL_AUTOCOMMIT)
    cursor = con.cursor()
    cursor.execute("CREATE TABLE IF NOT EXISTS Epitopes (epitopename VARCHAR(1000) NOT NULL, linearsequence VARCHAR(1000) NOT NULL, startingposition VARCHAR(1000), endingposition VARCHAR(1000), genebankid VARCHAR(1000), sourceorganism VARCHAR(1000), tcellid VARCHAR(1000), hostorganism VARCHAR(1000), eccelltype VARCHAR(1000), apccelltype VARCHAR(1000), mhc VARCHAR(1000), mhctype VARCHAR(1000), assaytype VARCHAR(1000), qualitativemeasurements VARCHAR(1000), quantitativemeasurements VARCHAR(1000));") #Create table if not exists
    f = open("database.sql", "rt", encoding="utf-8")
    y = 0
    while True:
        line = f.readline()
        if not line: break
        cursor.execute(line.strip("\n"))
        y += 1
        print(y)
    con.close()
    f.close()
    os.system("rm dumping.txt")


#Extract required information for training
def extract_from_psql(host1, database1, user1, password1, port1):
    #Extraact IFNg/IL-10 release
    con = psycopg2.connect(
                                host = host1,
                                database = database1,
                                user = user1 ,
                                password = password1,
                                port = port1
                                )
    con.set_isolation_level(ISOLATION_LEVEL_AUTOCOMMIT)
    cursor = con.cursor()
    #Obtain CSV
    fid = open(cytokine+'release.csv', 'w',encoding="utf-8")
    sql = "COPY (SELECT * FROM Epitopes WHERE mhctype = 'MHC II' AND assaytype = '"+cytokine+" release' and qualitativemeasurements = 'Positive') TO STDOUT WITH CSV HEADER;"
    cursor.copy_expert(sql, fid)
    con.close()

#Extract sequences and eliminate repited ones
def extract_sequences():
    f = open(cytokine+"release.csv", "rt", encoding="utf-8")
    list = []
    while True:
        line = f.readline()
        if not line: break
        sequence = line.split(",")[1]
        if sequence not in list:
            if sequence != "NA":
                if sequence != "linearsequence":
                    list.append(sequence)
    f.close()
    return list

def generate_training_file(list,seqNeg):
    f1 = open("training_"+cytokine+"release.txt", "wt", encoding="utf-8")
    f1.write("Peptide,Outcome")
    f1.write("\n")
    for x in list:
        f1.write(x+",1")
        f1.write("\n")
    for x in seqNeg:
        f1.write(x+",0")
        f1.write("\n")
    f1.close()

#Generate random binders
def generate_random_binders(pwd,list):
    os.system("mkdir CNN-PepPred/CNN-PepPred/Sequence")
    os.system("wget http://www.uniprot.org/uniprot/P98160.fasta")
    os.system("wget http://www.uniprot.org/uniprot/Q8WZ42.fasta")
    os.system("mv Q8WZ42.fasta P98160.fasta CNN-PepPred/CNN-PepPred/Sequence")
    os.system("mv CNN-PepPred/CNN-PepPred/* .")
    from generateRandomNonBinders import generateRandomNonBinders
    seqNeg = generateRandomNonBinders(pwd+"/Sequence", seq = list, N = len(list), maxFiles = 2)
    generate_training_file(list,seqNeg)

#Do the ML
def do_train_CV_logoplot(pwd,inpu1t):
    os.system("mkdir Results")
    f2 = open("train_CV_logoplot_apply.txt", "wt")
    f2.write("Input,Input_Value\nallele,"+cytokine+"release\nsavePath,"+pwd+"/Results\ndoTraining,1\ntrainingDataPath,"+pwd+"/training_"+cytokine+"release.txt\ndoLogoSeq,1\ndoCV,1\nkFold,\ndoApplyData,1\ntrainedModelsFile,\napplyDataPath,"+pwd+"/"+input1+"\nepitopesLength,15\nparametersFile,parameters.txt\nsaveClassObject,0")
    f2.close()
    import model_from_template
    modelCNN = model_from_template.main('train_CV_logoplot_apply.txt')

#Delete all non-necessary files
def delete_files():
    os.system("rm AssayTypeList.xml "+cytokine+"release.csv MhcAlleleNameList.xml OrganismList.xml blosum62.txt database.sql generateRandomNonBinders.py model_environment_cpu.yml model_environment_gpu.yml model_from_template.py model_initializer.py parameters.txt requirements_CPU.txt requirements_GPU.txt template1_Train_CV_logoPlot_Apply.txt template2_Train.txt template3_Apply.txt template_empty.txt template_pretrained_model_example.txt test_template.txt train_CV_logoplot_apply.txt training_"+cytokine+"release.txt")
    os.system("rm -rf CNN-PepPred Sequence __pycache__ trainedIEDBmodels Database Test trainedIEDBmodels_TL Example")

#Call functions
#Download files
importing_files()
#Dictionaries
org_dicc = obtain_dicc("OrganismList.xml", "OrganismId", "OrganismName")
assay_dicc = obtain_dicc("AssayTypeList.xml","AssayTypeId", "Response")
mhc_dicc = obtain_dicc("MhcAlleleNameList.xml", "MhcAlleleRestrictionId", "DisplayedRestriction")
#Iterate
iterate(assay_dicc, org_dicc ,mhc_dicc)
#Remove repeated lines
remove_lines()
#Enter PSQL, Create table, insert into table
enter_psql(host1, database1, user1, password1, port1)
#Obtain PWD
pwd = sp.getoutput('pwd')
#Extract required information from PSQL
extract_from_psql(host1, database1, user1, password1, port1)
#Extract sequences for training
list = extract_sequences()
#CNNPepPred
# Obtain CNNPepPred
os.system("git clone https://github.com/ComputBiol-IBB/CNN-PepPred")
#Generate random binders and generate training file
generate_random_binders(pwd, list)
#Do the machine learning
do_train_CV_logoplot(pwd,input1)
#Remove files
delete_files()

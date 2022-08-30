
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

#Read the input text
f = open("information.txt", "rt")
alelles = "NO"
pathogen = "NO"
while True:
    line = f.readline()
    if not line: break
    if "host:" in line:
        host1 = line.split(":")[1].strip("\n")
    if "database:" in line:
        database1 = line.split(":")[1].strip("\n")
    if "user:" in line:
        user1 = line.split(":")[1].strip("\n")
    if "password:" in line:
        password1 = line.split(":")[1].strip("\n")
    if "port:" in line:
        port1 = line.split(":")[1].strip("\n")
    if "input:" in line:
        if ".fasta" in line:
            input1 = line.split(":")[1].strip("\n")
        else:
            uniprot_ID = line.split(":")[1].strip("\n")
            print(uniprot_ID)
            os.system("wget http://www.uniprot.org/uniprot/"+uniprot_ID+".fasta")
            input1 = uniprot_ID+".fasta"
    if "prediction:" in line:
        cytokine_all = line.split(":")[1].strip("\n")
        cytokine_list = []
        for x in cytokine_all.split(","):
            cytokine_list.append(x)
    if "microorganism:" in line:
        if line.split(":")[1].strip("\n").lower() != "no":
            pathogen = line.split(":")[1].strip("\n")
    if "alelles:" in line:
        if line.split(":")[1].strip("\n").lower() == "yes":
            alelles = "YES"
    if "threshold:" in line:
        if line.split(":")[1].strip("\n").lower() != "no":
            threshold = float(line.split(":")[1].strip("\n"))
    if "length:" in line:
        length = str(line.split(":")[1].strip("\n"))
    if "update:" in line:
        update = str(line.split(":")[1].strip("\n")).lower()

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
                        line = "insert into Epitopes (epitopename, linearsequence, startingposition, endingposition, genebankid, sourceorganism, tcellid, hostorganism, eccelltype, apccelltype, mhc, mhctype, assaytype, qualitativemeasurements, quantitativemeasurements) values ("+"'"+ epitope+"' ,"+"'"+linearsequence + "' ,"+ "'"+sposition + "' ," + "'"+eposition+ "' ," + "'"+ genbankid + "' ," +"'"+ org_dicc[sourceorganism].replace("'","") + "' ," +"'"+cellid + "' ," +"'"+ org_dicc[hostorganism]  + "' ," + "'"+ ec_celltype+ "' ," + "'"+ apc_celltype+ "' ," + "'"+ mhc_dicc[mhc]+ "' ," + "'"+ mhc_class + "' ,"+ "'"+ assay_dicc[assaytype] +"' ,"  + "'"+ qualitativemeasurements+ "' ," + "'"+ quantitativemeasurements+"');"
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

#IF ALELLES
#Modify files
def modify_file(allele,directory):
    f = open("template3_Apply.txt","rt")
    f1 = open("inbetween.txt","w+")
    while True:
        line = f.readline()
        if not line: break
        match = re.search(r'allele,',line)
        match1 = re.search(r'trainedModelsFile',line)
        if match:
            f1.write(line[0:7]+allele+"\n")
        if match1:
            f1.write(line[0:18]+pwd+"/"+directory+"\n")
        if not match:
            if not match1:
                f1.write(line)
    f.close()
    f1.close()
    os.system("rm template3_Apply.txt")
    os.system("mv inbetween.txt template3_Apply.txt")

def obtain_alleles():
    # Obtain alleles
    for directory in list(glob.glob('trainedIEDBmodels/*')):
        allele = str(directory.strip("trainedIEDBmodels/model_"))
        print(directory)
        modify_file(allele,directory)
        import model_from_template
        modelCNN = model_from_template.main("template3_Apply.txt")

def obtain_peptidesfasta():
    #Obtain peptides.fasta
    alele = []
    list1 = []
    name = []
    beg = []
    end = []
    score = []
    for file in list(glob.glob("Alleles/*/*.txt")):
        os.system("cp "+file+" .")
    for file in list(glob.glob("HLA*")):
        y = 0
        f = open(file,"rt")
        while True:
            line = f.readline()
            if not line: break
            if line.split(",")[0] != "Peptide_Source":
                if float(line.split(",")[5].strip("\n")) >= threshold:
                    if line.split(",")[3] not in list1:
                        alele.append(file.strip("_predictedOutcome.txt"))
                        list1.append(line.split(",")[3])
                        name.append(line.split(",")[0])
                        beg.append(line.split(",")[1])
                        end.append(line.split(",")[2])
                        score.append(line.split(",")[5])
                        y += 1
        print(file.strip("_predictedOutcome.txt")+":",y," epitopes.")
        f.close()

    f1 = open("peptides.fasta","wt",encoding = "utf-8")
    for x in range(len(alele)):
        f1.write(">"+name[x]+">"+str(beg[x])+">"+str(end[x])+">"+str(score[x].strip("\n"))+">"+alele[x]+"\n")
        f1.write(list1[x]+"\n")
    f1.close()

#Connect to Postgress
def enter_psql(host1, database1, user1, password1, port1):
    #Connect to postgress
    con = psycopg2.connect(
                            host = host1,
                            database = database1,
                            user = user1 ,
                            password = password1,
                            port = port1
                            )
    con.set_isolation_level(ISOLATION_LEVEL_AUTOCOMMIT)
    cursor = con.cursor()
    #Create table if not exists
    cursor.execute("CREATE TABLE IF NOT EXISTS Epitopes (epitopename VARCHAR(10000), linearsequence VARCHAR(10000), startingposition VARCHAR(10000) , endingposition VARCHAR(10000) , genebankid VARCHAR(10000) , sourceorganism VARCHAR(10000) , tcellid VARCHAR(10000) , hostorganism VARCHAR(10000) , eccelltype VARCHAR(10000) , apccelltype VARCHAR(10000) , mhc VARCHAR(10000) , mhctype VARCHAR(10000) , assaytype VARCHAR(10000) , qualitativemeasurements VARCHAR(10000) , quantitativemeasurements VARCHAR(10000));")
    #Extract the current information of the table
    fid = open('whole_database.csv',  "wt", encoding = "utf-8")
    sql = "COPY (SELECT * FROM Epitopes) TO STDOUT WITH CSV HEADER;"
    cursor.copy_expert(sql, fid)
    fid.close()
    #List for appending all the information of the current table
    list = []
    y = 0#Counting
    #Open files
    f1 = open('whole_database.csv', "rt", encoding = "utf-8")
    f = open("database.sql", "rt", encoding="utf-8")
    #Append all rows in the list
    while True:
        line = f1.readline()
        if not line: break
        list.append(line.strip("\n"))
    f1.close()
    #Compare if the rows are already in the table
    while True:
        string = ""
        line = f.readline()
        if not line: break
        compare_line = line.strip("insert into Epitopes (epitopename, linearsequence, startingposition, endingposition, genebankid, sourceorganism, tcellid, hostorganism, eccelltype, apccelltype, mhc, mhctype, assaytype, qualitativemeasurements, quantitativemeasurements) values (")
        compare_line1 = compare_line.strip(");")
        for x in compare_line1.split(","):
            string += x[1:len(x)-2]+","
        compare_line2 = string[0:len(string)-1] #Obtain from the psql file a csv to compare with the list
        try:
            if compare_line2[0:len(compare_line2)-2].strip("\n") in list: #Comparing in the list
                print("Already in the table")
            else:
                y += 1
                print(y)
                cursor.execute(line.strip("\n"))
        except:
            pass
    f.close()
    con.close()
    os.system("rm whole_database.csv")
    os.system("rm dumping.txt")

#Extract required information for training
def extract_from_psql(host1, database1, user1, password1, port1, cytokine):
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
def extract_sequences(cytokine):
    f = open(cytokine.strip("/")+"release.csv", "rt", encoding="utf-8")
    list = []
    while True:
        line = f.readline()
        if not line: break
        sequence = line.split(",")[1]
        mic = line.split(",")[5]
        if sequence not in list:
            if sequence != "NA":
                if sequence != "linearsequence":
                    if pathogen not in mic:
                        list.append(sequence)
    f.close()
    return list

#Using epitopes from microorganism to validate
def epitopes_validating(cytokine):
    f = open(cytokine.strip("/")+'release.csv', 'rt',encoding="utf-8")
    f1 = open(cytokine.strip("/")+"epitopes.csv","wt", encoding = 'utf-8')
    list = []
    names = []
    while True:
        line = f.readline()
        if not line: break
        if pathogen in line:
            if cytokine in line:
                name = line.split(",")[0]
                sequence = line.split(",")[1]
                if sequence not in list:
                    if sequence != "NA":
                        if sequence != "linearsequence":
                            list.append(sequence)
                            names.append(name)
    for x in range(len(list)):
        f1.write(">"+names[x]+" \n")
        f1.write(list[x]+"\n")
    f.close()
    f1.close()
    #Unite all together
    os.system("cat "+cytokine.strip("/")+"epitopes.csv >> "+input1)
    os.system("rm "+cytokine.strip("/")+"epitopes.csv")
    #See if repeated:
    f2 = open(input1, 'rt',encoding="utf-8")
    input_epitopes = []
    input_names = []
    while True:
        line1 = f2.readline()
        line2 = f2.readline()
        if not line1: break
        if line2.strip("\n") not in input_epitopes:
            input_epitopes.append(line2.strip("\n"))
            input_names.append(line1.strip("\n"))
    f2.close()
    os.system("rm "+input1)
    #Rewrite the whole input without repetitions
    f2 = open(input1, 'wt',encoding="utf-8")
    for x in range(len(input_epitopes)):
        f2.write(input_names[x]+"\n")
        f2.write(input_epitopes[x]+"\n")
    f2.close()

#Generate training file
def generate_training_file(list,seqNeg,cytokine):
    f1 = open("training_"+cytokine.strip("/")+"release.txt", "wt", encoding="utf-8")
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
def generate_random_binders(pwd,list,cytokine):
    from generateRandomNonBinders import generateRandomNonBinders
    seqNeg = generateRandomNonBinders(pwd+"/Sequence", seq = list, N = len(list), maxFiles = 2)
    generate_training_file(list,seqNeg,cytokine)

#Do the ML with update
def do_train_CV_logoplot(pwd,input1,cytokine):
    f2 = open("train_CV_logoplot_apply.txt", "wt",encoding = "utf-8")
    f2.write("Input,Input_Value\nallele,"+cytokine+"release\nsavePath,"+pwd+"/Results\ndoTraining,1\ntrainingDataPath,"+pwd+"/training_"+cytokine.strip("/")+"release.txt\ndoLogoSeq,1\ndoCV,1\nkFold,\ndoApplyData,1\ntrainedModelsFile,\napplyDataPath,"+pwd+"/"+input1+"\nepitopesLength,"+length+"\nparametersFile,parameters.txt\nsaveClassObject,0")
    f2.close()
    import model_from_template
    modelCNN = model_from_template.main("train_CV_logoplot_apply.txt")

#Do the ML without update
def do_train_CV_logoplot_models(pwd,input1,cytokine):
    f2 = open("train_CV_logoplot_apply.txt", "wt",encoding = "utf-8")
    f2.write("Input,Input_Value\nallele,"+cytokine.replace("-","_")+"release\nsavePath,"+pwd+"/Results\ndoTraining,0\ntrainingDataPath,\ndoLogoSeq,0\ndoCV,0\nkFold,\ndoApplyData,1\ntrainedModelsFile,\napplyDataPath,"+pwd+"/"+input1+"\nepitopesLength,"+length+"\nparametersFile,parameters.txt\nsaveClassObject,0")
    f2.close()
    import model_from_template
    modelCNN = model_from_template.main("train_CV_logoplot_apply.txt")

#Delete all non-necessary files
def delete_files(cytokine):
    os.system("rm AssayTypeList.xml MhcAlleleNameList.xml OrganismList.xml blosum62.txt database.sql generateRandomNonBinders.py model_environment_cpu.yml model_environment_gpu.yml model_from_template.py model_initializer.py parameters.txt requirements_CPU.txt requirements_GPU.txt template1_Train_CV_logoPlot_Apply.txt template2_Train.txt template3_Apply.txt template_empty.txt template_pretrained_model_example.txt test_template.txt train_CV_logoplot_apply.txt")
    os.system("rm -rf CNN-PepPred Sequence __pycache__ trainedIEDBmodels Database Test trainedIEDBmodels_TL Example")


#Call functions
#Download files
if update == "yes":
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
    #Download all required files
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

#Get PWD
pwd = sp.getoutput('pwd')

#Obtain resulsts for all human aleles using the models
if alelles == "YES":
    os.system("git clone https://github.com/ComputBiol-IBB/CNN-PepPred")
    os.system("mv CNN-PepPred/CNN-PepPred/* .")
    os.system("mkdir Alleles")
    os.system("rm template3_Apply.txt")
    template = open("template3_Apply.txt","wt",encoding = "utf-8")
    template.write("Input,Input_Value\nallele,HLA_DQA1_01_02_DQB1_06_02\nsavePath,"+pwd+"/Alleles\ndoTraining,0\ntrainingDataPath,\ndoLogoSeq,0\ndoCV,0\nkFold,\ndoApplyData,1\ntrainedModelsFile,\napplyDataPath,"+pwd+"/"+input1+"\nepitopesLength,"+length+"\nparametersFile,parameters.txt\nsaveClassObject,0")
    template.close()

    #Calling functions
    obtain_alleles()
    obtain_peptidesfasta()
    os.system("rm blosum62.txt generateRandomNonBinders.py model_environment_cpu.yml model_environment_gpu.yml model_from_template.py model_initializer.py parameters.txt requirements_CPU.txt requirements_GPU.txt template1_Train_CV_logoPlot_Apply.txt template2_Train.txt template3_Apply.txt template_empty.txt template_pretrained_model_example.txt test_template.txt train_CV_logoplot_apply.txt training_")
    os.system("rm -rf CNN-PepPred Sequence __pycache__ trainedIEDBmodels Database Test trainedIEDBmodels_TL Example")
    input1 = "peptides.fasta"
    os.system("rm *_predictedOutcome.txt")

# Obtain CNNPepPred and prepare for running
os.system("git clone https://github.com/ComputBiol-IBB/CNN-PepPred")
os.system("mkdir CNN-PepPred/CNN-PepPred/Sequence")
os.system("wget http://www.uniprot.org/uniprot/P98160.fasta")
os.system("wget http://www.uniprot.org/uniprot/Q8WZ42.fasta")
os.system("mv Q8WZ42.fasta P98160.fasta CNN-PepPred/CNN-PepPred/Sequence")
os.system("mv CNN-PepPred/CNN-PepPred/* .")
os.system("mkdir Results")

#Prediction for cytokines
for x in cytokine_list:
    if update == "yes":
        #Extract required information from PSQL
        extract_from_psql(host1, database1, user1, password1, port1, x)
        #Extract sequences for training
        list = extract_sequences(x)
        #Insert epitopes to input for validation of the programme
        if pathogen != "NO":
            epitopes_validating(x)
        #Generate random binders and generate training file
        generate_random_binders(pwd, list, x)
        #Do the machine learning
        print("TRAINING OF: "+x)
        do_train_CV_logoplot(pwd,input1,x)
        os.system("rm "+x.strip("/")+"release.csv training_"+x.strip("/")+"release.txt")
    if update == "no":
        os.system("cp "+pwd+"/Cytokines_models/* "+pwd+"/trainedIEDBmodels")
        #Do the machine learning
        print("TRAINING OF: "+x)
        do_train_CV_logoplot_models(pwd,input1,x)

#Remove files all not-wanted files
delete_files(x)

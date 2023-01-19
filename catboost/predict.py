import argparse
import json
import catboost
from catboost import CatBoostClassifier
import ROOT
import numpy as np
import uproot3 as uproot
import awkward as ak
import os 
import array
import yaml

parser = argparse.ArgumentParser()

parser.add_argument('--raw_data_dir' , type=str, default = "/volatile/clas12/users/gmat/clas12analysis.sidis.data/rga/ML/raw", help="Location of raw ROOT files for predicting [default: /volatile/clas12/users/gmat/clas12analysis.sidis.data/rga/ML/raw]")
parser.add_argument('--preprocess_data_dir' , type=str, default = "/volatile/clas12/users/gmat/clas12analysis.sidis.data/rga/ML/catboost/preprocess_pi0", help="Location of preprocessed ROOT files for predicting [default: /volatile/clas12/users/gmat/clas12analysis.sidis.data/rga/ML/catboost/preprocess_pi0]")
parser.add_argument('--output_dir' , type=str, default= "/volatile/clas12/users/gmat/clas12analysis.sidis.data/rga/ML/catboost/predict_pi0", help="Location of output ROOT TFiles containing predictions [default: /volatile/clas12/users/gmat/clas12analysis.sidis.data/rga/ML/catboost/predict_pi0]")
parser.add_argument('--model_dir' , type=str, default='/work/clas12/users/gmat/nnPi0/catboost/models/model', help='full path to directory for model and plots [default: /work/clas12/users/gmat/nnPi0/catboost/models/model]')
parser.add_argument('--subdata' , type=str, default="all" ,help="Specifies the subset of data files to use for predicting [default: all] (SEE ./utils/subdata.json FOR OPTIONS)")
parser.add_argument('--version' , type=str, default="MC" , help="Type of dataset (MC or nSidis) [default: MC]")

FLAGS = parser.parse_args()
RAW_DATA_DIR    = FLAGS.raw_data_dir
PREPROCESS_DATA_DIR    = FLAGS.preprocess_data_dir
OUTPUT_DIR  = FLAGS.output_dir
MODEL_DIR  = FLAGS.model_dir
VERSION     = FLAGS.version
SUBDATA    = FLAGS.subdata
INPUT_YAML = MODEL_DIR+"/inputs.yaml"
with open(INPUT_YAML,'r') as file:
    inyaml = yaml.safe_load(file)

BRANCH_NAMES = ['flag','ievent']  +  inyaml["inputs"]

fjs = open ('utils/subdata.json', "r")
JSON = json.loads(fjs.read())
SUBDATA_KEYS=[key for key in JSON.keys()]
assert(SUBDATA=="all" or SUBDATA in SUBDATA_KEYS)

MODEL_FILE=MODEL_DIR+"/catboost_model"

RAW_FILES = [] # Contains event-by-event information
PREPROCESSED_FILES = [] # Contains photon-by-photon information
POSTPROCESSED_FILES = [] # Contains event-by-event information with CatBoost weights

assert(VERSION=="nSidis" or VERSION=="MC")

def load_model():
    from_file = CatBoostClassifier()
    from_file.load_model(MODEL_FILE)
    return from_file

def load_files():

    for file in os.listdir(PREPROCESS_DATA_DIR):
        if (file.endswith(".root") and VERSION in file):
            foundFile=False
            if(SUBDATA!="all"):
                for RUN in JSON[SUBDATA]:
                    if RUN in file:
                        foundFile=True
            else:
                foundFile=True
            if(not foundFile):
                continue
            if(not "5032" in file):
                continue
            RAW_FILES.append(RAW_DATA_DIR+"/"+file)
            PREPROCESSED_FILES.append(PREPROCESS_DATA_DIR+"/"+file+":PreProcessedEvents")
            POSTPROCESSED_FILES.append(OUTPUT_DIR+"/"+file)
            
    print("\n",len(RAW_FILES),"root files found for the ML predictions\n","-"*70,"\n")

def load_data(file):

    #load the tree
    tree = uproot.open(file.split(":")[0]) # Because of uproot4 shenanigans
    tree = tree["PreProcessedEvents"]
    
    # Column idx
    idx_ievent = BRANCH_NAMES.index("ievent")
    idx_flag   = BRANCH_NAMES.index("flag")

    #load the branches into a numpy array
    data = np.array([tree[b].array() for b in BRANCH_NAMES], dtype=np.float32).T

    #reshape the data into a N by (number of TBranches) matrix
    data = data.reshape(data.shape[0],len(BRANCH_NAMES))
    
    # Delete ievent and flag columns
    X =  np.delete(data,[idx_ievent,idx_flag],1)
    y =  data[:,idx_flag]
    iev = data[:,idx_ievent].astype(int)
    
    return X,y,np.unique(iev)

def predict():
    
    load_files()
    
    model = load_model()
    
    ifile=0
    for rawfile,prefile,postfile in zip(RAW_FILES,PREPROCESSED_FILES,POSTPROCESSED_FILES):
        ifile+=1
        print("Predicting file",ifile,"of",len(RAW_FILES))
        print("Loading data...",end="")
        X,y,ievent = load_data(prefile)
        print("Complete")
        print("Making predictions...",end="")
        #X["prob"]=model.predict_proba(X)[:,1]
        prob=model.predict_proba(X)[:,1]
        print("Complete")
        # Open the existing TTree
        raw_tfile = ROOT.TFile(rawfile, "READ")
        raw_tree = raw_tfile.Get("RawEvents")
#         u=uproot.open(rawfile)
#         u=u["RawEvents"]
#         pid=u.arrays("pid")
#         print(pid[0])
        #pid=ak.to_list(pid)
        
        #pid=u["pid"].array(library="np")
        post_tfile = ROOT.TFile(postfile,"RECREATE")
        post_tree = raw_tree.CloneTree(0)
        post_tree.SetName("PostProcessedEvents")
        listOfBranches = post_tree.GetListOfBranches()
        
        #Get the list of branches
        branchList = post_tree.GetListOfBranches()

        #Loop over the list of branches
        for branch in branchList:
            #Check if the branch is in the list
            if branch.GetName() not in listOfBranches:
                #Remove the branch from the tree
                post_tree.SetBranchStatus(branch.GetName(), 0)

        #create a new branch for the tree
        weights = array.array('f',100*[0.0])
        post_tree.Branch("catboost_weight", weights,'catboost_weight[nPart]/F')
        
        # Fill the new branch with data
        print("Creating PostProcess Tree...",end="")
        k=0
        
        for i in ievent:
            raw_tree.GetEntry(i)
            pid=np.array(raw_tree.pid)
            for j,PID in enumerate(pid):
                if(PID==22):
                    weights[j]=prob[k]
                    k+=1
                else:
                    weights[j]=1
            post_tree.Fill()
        # Write the updated tree to a new file
        print("Complete\n")
        post_tree.Write()
        post_tfile.Close()
        raw_tfile.Close()
    
    print("\n All done \n")
    
if __name__ == "__main__":
    predict()

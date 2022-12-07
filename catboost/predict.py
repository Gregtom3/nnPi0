import argparse
import catboost
from catboost import CatBoostClassifier
import ROOT
import numpy as np
import uproot
import os 
import array

parser = argparse.ArgumentParser()

parser.add_argument('--raw_data_dir' , type=str, default = "/volatile/clas12/users/gmat/clas12analysis.sidis.data/rga/ML/raw", help="Location of raw ROOT files for predicting [default: /volatile/clas12/users/gmat/clas12analysis.sidis.data/rga/ML/raw]")
parser.add_argument('--preprocess_data_dir' , type=str, default = "/volatile/clas12/users/gmat/clas12analysis.sidis.data/rga/ML/preprocess_catboost", help="Location of preprocessed ROOT files for predicting [default: /volatile/clas12/users/gmat/clas12analysis.sidis.data/rga/ML/preprocess_catboost]")
parser.add_argument('--output_dir' , type=str, default= "/volatile/clas12/users/gmat/clas12analysis.sidis.data/rga/ML/postprocess_catboost", help="Location of output ROOT TFiles containing predictions [default: /volatile/clas12/users/gmat/clas12analysis.sidis.data/rga/ML/postprocess_catboost]")
parser.add_argument('--model_file', type=str, default='./models/test', help='Path to trained CatBoost model [default: ./models/test]')
parser.add_argument('--version' , type=str, default="MC" , help="Type of dataset (MC or nSidis) [default: MC]")

FLAGS = parser.parse_args()
RAW_DATA_DIR    = FLAGS.raw_data_dir
PREPROCESS_DATA_DIR    = FLAGS.preprocess_data_dir
OUTPUT_DIR  = FLAGS.output_dir
MODEL_FILE  = FLAGS.model_file
VERSION     = FLAGS.version

RAW_FILES = [] # Contains event-by-event information
PREPROCESSED_FILES = [] # Contains photon-by-photon information
POSTPROCESSED_FILES = [] # Contains event-by-event information with CatBoost weights

assert(VERSION=="nSidis" or VERSION=="MC")

def load_model():
    from_file = CatBoostClassifier()
    from_file.load_model(MODEL_FILE)
    return from_file

def load_files():

    for file in os.listdir(RAW_DATA_DIR):
        if (file.endswith(".root") and VERSION in file):
            RAW_FILES.append(RAW_DATA_DIR+"/"+file)
            PREPROCESSED_FILES.append(PREPROCESS_DATA_DIR+"/"+file+":PreProcessedEvents")
            POSTPROCESSED_FILES.append(OUTPUT_DIR+"/"+file)
    print(len(RAW_FILES),"root files found for the ML predictions")
    
def load_data(file):

    #load the tree
    tree = uproot.open(file)

    #get the branch names
    branch_names = tree.keys()

    # Column idx
    idx_ievent = branch_names.index("ievent")
    idx_flag   = branch_names.index("flag")

    #load the branches into a numpy array
    data = np.array([tree[b].array() for b in branch_names], dtype=np.float32).T

    #reshape the data into a N by (number of TBranches) matrix
    data = data.reshape(data.shape[0],len(branch_names))
    
    # Delete ievent and flag columns
    X =  np.delete(data,[idx_ievent,idx_flag],1)
    y =  data[:,idx_flag]
    iev =data[:,idx_ievent]
    
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

        post_tfile = ROOT.TFile(postfile,"RECREATE")
        post_tree = raw_tree.CloneTree(0)
        post_tree.SetName("PostProcessedEvents")

        #listOfBranches is the list of desired branches
        listOfBranches = ["hel","x","Q2","W","nPart","px","py","pz","E","pid","theta","eta","phi"]

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
        for i in range(len(ievent)):
            raw_tree.GetEntry(i)

            nPart   = raw_tree.nPart
            
            for j in range(nPart):
                pid = raw_tree.pid[j]
                if(pid==22):
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

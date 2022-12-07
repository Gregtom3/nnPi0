import argparse
import ROOT
import pandas as pd
import numpy as np
import uproot
import os

from sklearn.model_selection import train_test_split
from catboost import CatBoostClassifier, Pool, metrics, cv
from sklearn.metrics import accuracy_score
import catboost
from catboost.utils import get_confusion_matrix

parser = argparse.ArgumentParser()

parser.add_argument('--lr' , type=float, default = "0.05", help="Learning rate [default: 0.05]")
parser.add_argument('--depth', type=int, default = "5", help="Maximum depth per tree [default: 5]")
parser.add_argument('--l2_leaf_reg', type=float, default = "4.0", help="L2 leaf reg [default: 4.0]")
parser.add_argument('--iters', type=int, default="1000", help="Number of iterations [default: 1000]")
parser.add_argument('--max_leaves', type=int, default="64", help="Maximum number of leaves [default: 64]")
parser.add_argument('--min_data_in_leaf', type=int, default="1", help="Minimum number of data for leaf to split [default: 1]")
parser.add_argument('--train_size', type=float, default="0.75", help="Fraction of sample to use for training [default: 0.75]")
parser.add_argument('--data_dir', type=str, default='./', help='directory with data [default: ./]')
parser.add_argument('--seed' , type=int , default="42" , help="Random seed [default: 42]")

FLAGS = parser.parse_args()
LR    = FLAGS.lr
DEPTH = FLAGS.depth
L2_LEAF_REG = FLAGS.l2_leaf_reg
ITERS = FLAGS.iters
MAX_LEAVES = FLAGS.max_leaves
MIN_DATA_IN_LEAF = FLAGS.min_data_in_leaf
TRAIN_SIZE = FLAGS.train_size
DATA_DIR = FLAGS.data_dir
SEED  = FLAGS.seed

if not os.path.exists(DATA_DIR):
    print("ERROR:",DATA_DIR,"not found...Aborting...")

def get_data():

    root_files = []

    for file in os.listdir(DATA_DIR):
        if file.endswith(".root"):
            root_files.append(DATA_DIR+"/"+file+":PreProcessedEvents")

    print(len(root_files),"root files found for the ML train/test")

    #load the tree
    tree = uproot.open(root_files[0])

    #get the branch names
    branch_names = ['flag','ievent','nPhotons',
            'nHadrons',
             'gE',
             'gTheta',
             'gPhi',
             'g_pcal_e',
             'g1_pcal_e',
             'g2_pcal_e',
             'g_pcal_du',
             'g_pcal_dv',
             'g_pcal_m2u',
             'g_pcal_m2v',
             'g_pcal_m3u',
             'g_pcal_m3v',
             'g1R',
             'g2R',
             'g1M',
             'g2M',
             'g1dE',
             'g2dE',
             'h1R',
             'h2R',
             'h1M',
             'h2M',
             'h1dE',
             'h2dE',
             'h1q',
             'h2q',
             'eR',
             'eM',
             'edE']
    
    # Column idx
    idx_ievent = branch_names.index("ievent")
    idx_flag   = branch_names.index("flag")
    
    data = np.empty((0, len(branch_names)))
    
    # loop over the input tfiles
    for tfile in root_files:
    
        # open uproot TTree
        tree = uproot.open(tfile)

        #load the branches into a numpy array
        temp_data = np.array([tree[b].array() for b in branch_names], dtype=np.float32).T

        # add the numpy array to the overall array
        data = np.vstack([data, temp_data])
        break
    #reshape the data into a N by (number of TBranches) matrix
    data = data.reshape(data.shape[0],len(branch_names))
    
    # Delete ievent and flag columns
    X = np.delete(data,[idx_ievent,idx_flag],1)
    y=  data[:,idx_flag]
    
    X_train, X_validation, y_train, y_validation = train_test_split(X, y, train_size=TRAIN_SIZE, random_state=SEED)
    return X_train, X_validation, y_train, y_validation

def train():
    
    X_train, X_validation, y_train, y_validation = get_data()
    
    numeric_train_pool = Pool(X_train, y_train)
    numeric_val_pool = Pool(X_validation, y_validation)
    
    ops = {'learning_rate': LR,
               'depth': DEPTH,
               'l2_leaf_reg': L2_LEAF_REG,
               'random_seed': SEED,
               'iterations': ITERS,
               #'max_leaves': MAX_LEAVES,
               'min_data_in_leaf': MIN_DATA_IN_LEAF
        }

    model = CatBoostClassifier(**ops,
                            custom_loss=[metrics.Accuracy()], 
                            task_type="GPU",
                            devices='0:1')
    model.fit(numeric_train_pool, verbose=1, eval_set=numeric_val_pool)
    
    model.save_model("models/test")
    
    cm = get_confusion_matrix(model, numeric_val_pool)
    print(cm)
    
if __name__ == "__main__":
    train()

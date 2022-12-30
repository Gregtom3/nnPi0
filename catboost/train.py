import argparse
import ROOT
import pandas as pd
import numpy as np
import uproot
import os
import matplotlib.pyplot as plt
import seaborn as sns
import shutil
import json
import yaml

from sklearn.model_selection import train_test_split
from catboost import CatBoostClassifier, Pool, metrics, cv
from sklearn.metrics import accuracy_score
import catboost
from catboost.utils import get_confusion_matrix
from catboost.utils import get_roc_curve, select_threshold, get_fpr_curve

plt.style.use('science')
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})
parser = argparse.ArgumentParser()

parser.add_argument('--lr' , type=float, default = "0.05", help="Learning rate [default: 0.05]")
parser.add_argument('--depth', type=int, default = "5", help="Maximum depth per tree [default: 5]")
parser.add_argument('--l2_leaf_reg', type=float, default = "4.0", help="L2 leaf reg [default: 4.0]")
parser.add_argument('--iters', type=int, default="1000", help="Number of iterations [default: 1000]")
parser.add_argument('--max_leaves', type=int, default="64", help="Maximum number of leaves [default: 64]")
parser.add_argument('--min_data_in_leaf', type=int, default="1", help="Minimum number of data for leaf to split [default: 1]")
parser.add_argument('--train_size', type=float, default="0.75", help="Fraction of sample to use for training [default: 0.75]")
parser.add_argument('--data_dir', type=str, default='./', help='directory with data [default: ./]')
parser.add_argument('--subdata', type=str, default='all', help='Specifies the MC files to be used in training (for specifying inbending vs outbending sets) [default: all] (SEE ./utils/subdata.json FOR OPTIONS)')
parser.add_argument('--model_dir' , type=str, default='/work/clas12/users/gmat/nnPi0/catboost/models/model', help='full path to directory for model and plots [default: /work/clas12/users/gmat/nnPi0/catboost/models/model]')
parser.add_argument('--input_yaml' , type=str, default="/work/clas12/users/gmat/nnPi0/catboost/input/input_noresonance.yaml",help="YAML file for model inputs [default: /work/clas12/users/gmat/nnPi0/catboost/input/input_noresonance.yaml]")
parser.add_argument('--make_plots', type=bool, default='false', help='create model performance plots in model_dir [default: false]')
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
SUBDATA  = FLAGS.subdata
SEED  = FLAGS.seed
MODEL_DIR = FLAGS.model_dir
MAKE_PLOTS = FLAGS.make_plots
INPUT_YAML = FLAGS.input_yaml
with open(INPUT_YAML,'r') as file:
    inyaml = yaml.safe_load(file)

BRANCH_NAMES = ['flag','ievent']  +  inyaml["inputs"]
FEATURE_LIST=BRANCH_NAMES[2:]


fjs = open ('utils/subdata.json', "r")
JSON = json.loads(fjs.read())
SUBDATA_KEYS=[key for key in JSON.keys()]
assert(SUBDATA=="all" or SUBDATA in SUBDATA_KEYS)

if os.path.exists(MODEL_DIR):
    shutil.rmtree(MODEL_DIR)
os.mkdir(MODEL_DIR)

# Copy the YAML file into the MODEL DIR
shutil.copyfile(INPUT_YAML,MODEL_DIR+"/inputs.yaml")

if not os.path.exists(DATA_DIR):
    print("ERROR:",DATA_DIR,"not found...Aborting...")

def get_data():

    root_files = []

    for file in os.listdir(DATA_DIR):
        if (file.endswith(".root") and "MC" in file):
            foundFile=False
            if(SUBDATA!="all"):
                for RUN in JSON[SUBDATA]:
                    if RUN in file:
                        foundFile=True
            else:
                foundFile=True
            if(foundFile):
                root_files.append(DATA_DIR+"/"+file+":PreProcessedEvents")

    print(len(root_files),"root files found for the ML train/test")

    #load the tree
    tree = uproot.open(root_files[0])

    #get the branch names
    branch_names = BRANCH_NAMES
    
    # Column idx
    idx_ievent = branch_names.index("ievent")
    idx_flag   = branch_names.index("flag")
    
    data = np.empty((0, len(branch_names)))
    
    FEATURE_LIST = branch_names[2:]
    # loop over the input tfiles
    for tfile in root_files:
    
        # open uproot TTree
        tree = uproot.open(tfile)

        #load the branches into a numpy array
        temp_data = np.array([tree[b].array() for b in branch_names], dtype=np.float32).T

        # add the numpy array to the overall array
        data = np.vstack([data, temp_data])

    #reshape the data into a N by (number of TBranches) matrix
    data = data.reshape(data.shape[0],len(branch_names))
    
    # Delete ievent and flag columns
    X = np.delete(data,[idx_ievent,idx_flag],1)
    y=  data[:,idx_flag]
    
    X_train, X_validation, y_train, y_validation = train_test_split(X, y, train_size=TRAIN_SIZE, random_state=SEED)
    return X_train, X_validation, y_train, y_validation


# SEE https://medium.com/@dtuk81/confusion-matrix-visualization-fc31e3f30fea for reference
def make_confusion_matrix(cf,
                          group_names=None,
                          categories='auto',
                          count=True,
                          percent=True,
                          cbar=True,
                          xyticks=True,
                          xyplotlabels=True,
                          sum_stats=True,
                          figsize=None,
                          cmap='Blues',
                          title=None):
    '''
    This function will make a pretty plot of an sklearn Confusion Matrix cm using a Seaborn heatmap visualization.
    Arguments
    ---------
    cf:            confusion matrix to be passed in
    group_names:   List of strings that represent the labels row by row to be shown in each square.
    categories:    List of strings containing the categories to be displayed on the x,y axis. Default is 'auto'
    count:         If True, show the raw number in the confusion matrix. Default is True.
    normalize:     If True, show the proportions for each category. Default is True.
    cbar:          If True, show the color bar. The cbar values are based off the values in the confusion matrix.
                   Default is True.
    xyticks:       If True, show x and y ticks. Default is True.
    xyplotlabels:  If True, show 'True Label' and 'Predicted Label' on the figure. Default is True.
    sum_stats:     If True, display summary statistics below the figure. Default is True.
    figsize:       Tuple representing the figure size. Default will be the matplotlib rcParams value.
    cmap:          Colormap of the values displayed from matplotlib.pyplot.cm. Default is 'Blues'
                   See http://matplotlib.org/examples/color/colormaps_reference.html
                   
    title:         Title for the heatmap. Default is None.
    '''


    # CODE TO GENERATE TEXT INSIDE EACH SQUARE
    blanks = ['' for i in range(cf.size)]

    if group_names and len(group_names)==cf.size:
        group_labels = ["{}\n".format(value) for value in group_names]
    else:
        group_labels = blanks

    if count:
        group_counts = ["{0:0.0f}\n".format(value) for value in cf.flatten()]
    else:
        group_counts = blanks

    if percent:
        group_percentages = ["{0:.2%}".format(value) for value in cf.flatten()/np.sum(cf)]
    else:
        group_percentages = blanks

    box_labels = [f"{v1}{v2}{v3}".strip() for v1, v2, v3 in zip(group_labels,group_counts,group_percentages)]
    box_labels = np.asarray(box_labels).reshape(cf.shape[0],cf.shape[1])


    # CODE TO GENERATE SUMMARY STATISTICS & TEXT FOR SUMMARY STATS
    if sum_stats:
        #Accuracy is sum of diagonal divided by total observations
        accuracy  = np.trace(cf) / float(np.sum(cf))

        #if it is a binary confusion matrix, show some more stats
        if len(cf)==2:
            #Metrics for Binary Confusion Matrices
            precision = cf[1,1] / sum(cf[:,1])
            recall    = cf[1,1] / sum(cf[1,:])
            f1_score  = 2*precision*recall / (precision + recall)
            stats_text = "\n\nAccuracy={:0.3f}\nPrecision={:0.3f}\nRecall={:0.3f}\nF1 Score={:0.3f}".format(
                accuracy,precision,recall,f1_score)
        else:
            stats_text = "\n\nAccuracy={:0.3f}".format(accuracy)
    else:
        stats_text = ""


    # SET FIGURE PARAMETERS ACCORDING TO OTHER ARGUMENTS
    if figsize==None:
        #Get default figure size if not set
        figsize = plt.rcParams.get('figure.figsize')

    if xyticks==False:
        #Do not show categories if xyticks is False
        categories=False


    # MAKE THE HEATMAP VISUALIZATION
    plt.figure(figsize=figsize,dpi=150)
    sns.heatmap(cf,annot=box_labels,fmt="",cmap=cmap,cbar=cbar,xticklabels=categories,yticklabels=categories)

    if xyplotlabels:
        plt.ylabel('True label')
        plt.xlabel('Predicted label' + stats_text)
    else:
        plt.xlabel(stats_text)
    
    if title:
        plt.title(title)
    
    plt.savefig(MODEL_DIR+"/confusion_matrix.pdf")
    
def make_fpr_curve(numeric_val_pool=0,model=0,FPR=0.01,myLabel=""):
    
    appTitle=""
    if(myLabel!=""):
        appTitle="({})".format(myLabel)
    roc_curve_values=get_roc_curve(model,numeric_val_pool)  
    x,y=get_fpr_curve(curve=roc_curve_values)
    plt.figure(dpi=150)
    plt.plot(x,y,color="black")
    plt.grid()
    plt.xlabel("Threshold (p)")
    plt.ylabel("False Positive Rate")
    plt.title("FPR Curve ".format(appTitle))
    plt.xlim(0,1)

    boundary = select_threshold(model, 
                                curve=roc_curve_values,  
                                FPR=FPR)
    if(appTitle):
        plt.text(0.2,0.64,"Threshold for FPR = {}\% : {:.3f}".format(FPR*100,boundary),fontsize=9)
        plt.yscale("log")
        plt.plot([0,boundary],[FPR,FPR],color="red")
        plt.plot([boundary,boundary],[y[1],FPR],color="red")
        plt.savefig(MODEL_DIR+"/fpr_curve_{}.pdf".format(myLabel))
    else:
        plt.savefig(MODEL_DIR+"/fpr_curve.pdf")
    
    roc_curve_values=get_roc_curve(model,numeric_val_pool)
    x=roc_curve_values[0]
    y=roc_curve_values[1]
    TPR=y[np.abs(x-FPR).argmin()]
    plt.figure(dpi=150)
    plt.plot(x,y,"k")
    plt.grid()
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("ROC Curve {}".format(appTitle))
    if(appTitle):
        plt.text(0.2,0.64,"TPR({}) = {}\% ".format(FPR,np.round(TPR*100,2)),fontsize=9)
        plt.plot([1,FPR],[TPR,TPR],color="red")
        plt.plot([FPR,FPR],[TPR,0],color="red")
        plt.savefig(MODEL_DIR+"/roc_curve_{}.pdf".format(myLabel))
    else:
        plt.savefig(MODEL_DIR+"/roc_curve.pdf".format(myLabel))

def save_param_importance(model):
    dfPars = pd.DataFrame(data={"Parameter": FEATURE_LIST,"Importance": model.get_feature_importance()})
    dfPars = dfPars.sort_values(by="Importance",ascending=False)
    dfPars.to_csv(MODEL_DIR+"/param_importance.csv",index=False)
    
def savefigs(numeric_val_pool, model):
    
    cm = get_confusion_matrix(model, numeric_val_pool)
    categories = ['Bkg','Signal']
    make_confusion_matrix(cm, 
                          categories=categories,
                          cmap='Blues',
                         cbar=False,
                         figsize=(4,3),
                         count=False,
                         title="Confusion Matrix\n on validation set")
    
    make_fpr_curve(numeric_val_pool,model,0.01,"Strict")
    make_fpr_curve(numeric_val_pool,model,0.03,"Medium")
    make_fpr_curve(numeric_val_pool,model,0.1,"Loose")
    make_fpr_curve(numeric_val_pool,model)
    
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
    
    model.save_model(MODEL_DIR+"/catboost_model")
    
    if(MAKE_PLOTS):
        savefigs(numeric_val_pool,model)
    
    save_param_importance(model)
        
if __name__ == "__main__":
    train()

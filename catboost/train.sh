python train.py --data_dir /volatile/clas12/users/gmat/clas12analysis.sidis.data/rga/ML/catboost/preprocess_pi0 \
    --depth 12 \
    --make_plots true \
    --subdata MC_inbending \
    --model_dir inbending 
    
python train.py --data_dir /volatile/clas12/users/gmat/clas12analysis.sidis.data/rga/ML/catboost/preprocess_pi0 \
    --depth 12 \
    --make_plots true \
    --subdata MC_outbending \
    --model_dir outbending 
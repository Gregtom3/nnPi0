python train.py --data_dir /volatile/clas12/users/gmat/clas12analysis.sidis.data/rga/ML/catboost/preprocess_pi0 \
    --depth 12 \
    --make_plots true \
    --subdata MC_inbending \
    --model_dir inbending_noresonance \
    --input_yaml ./input/input_noresonance.yaml
    
python train.py --data_dir /volatile/clas12/users/gmat/clas12analysis.sidis.data/rga/ML/catboost/preprocess_pi0 \
    --depth 12 \
    --make_plots true \
    --subdata MC_outbending \
    --model_dir outbending_noresonance \
    --input_yaml ./input/input_noresonance.yaml
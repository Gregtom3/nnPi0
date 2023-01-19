python train.py --data_dir /volatile/clas12/users/gmat/clas12analysis.sidis.data/rga/ML/catboost/preprocess_pi0 \
    --depth 12 \
    --make_plots true \
    --subdata MC_inbending \
    --model_dir /work/clas12/users/gmat/nnPi0/catboost/models/inbending/ \
    --input_yaml /work/clas12/users/gmat/nnPi0/catboost/input/input_noresonance.yaml
    
python train.py --data_dir /volatile/clas12/users/gmat/clas12analysis.sidis.data/rga/ML/catboost/preprocess_pi0 \
    --depth 12 \
    --make_plots true \
    --subdata MC_outbending \
    --model_dir /work/clas12/users/gmat/nnPi0/catboost/models/outbending \
    --input_yaml /work/clas12/users/gmat/nnPi0/catboost/input/input_noresonance.yaml
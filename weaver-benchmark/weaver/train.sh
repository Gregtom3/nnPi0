PREFIX='particlenet-dihadron-simple'
MODEL_CONFIG='particlenet_6pf.py'
DATA_CONFIG='simpleDihadron_pf_points_features.yaml'
PATH_TO_SAMPLES='/volatile/clas12/users/gmat/clas12analysis.sidis.data/rga/ML'

python train.py \
 --data-train ${PATH_TO_SAMPLES}'/tinyTrain.root' \
 --data-val ${PATH_TO_SAMPLES}'/tinyTest.root' \
 --fetch-by-file --fetch-step 1 --num-workers 3 \
 --data-config 'data/'${DATA_CONFIG} \
 --network-config networks/${MODEL_CONFIG} \
 --model-prefix output/${PREFIX} \
 --gpus 0,1 --batch-size 1024 --start-lr 5e-3 --num-epochs 10 --optimizer ranger \
 --log output/${PREFIX}.train.log

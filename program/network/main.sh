#!/bin/bash

current_path=$(pwd)
echo ${current_path}
cd ..
cd ..
main_path=$(pwd)
cd ${current_path}
echo $(pwd)

train_path=${main_path}/data/dataset/cross_val
test_file=${main_path}/data/dataset/independent_test/independent_test.csv
result_path=${main_path}/data/result

kfold=5
seqwin=41
machine_method=LGBM # RF SVM XGB LGBM

for machine_method in LGBM 
do
for encode_method in RCKmer DNC TNC CKSNAP ANF PseEIIP # ENAC binary NCP EIIP RCKmer DNC TNC CKSNAP ANF PseEIIP
do 
echo ${machine_method}
echo ${encode_method}

python ml_train_test_2OM.py  --intrain ${train_path} --intest ${test_file} --outpath ${result_path} --machine ${machine_method}  --encode ${encode_method} --fold ${kfold} --seqwin ${seqwin}

done
done


#!/bin/bash

ERA=$1

# Samples Run2016
ARTUS_OUTPUTS_2016=/ceph/swozniewski/SM_Htautau/ntuples/Artus16_Nov/Artus_2016_Jan19/merged_fixedjets2
ARTUS_FRIENDS_ET_2016=/ceph/swozniewski/SM_Htautau/ntuples/Artus16_Nov/Artus_2016_Jan19/et_keras_20190124
ARTUS_FRIENDS_MT_2016=/ceph/swozniewski/SM_Htautau/ntuples/Artus16_Nov/Artus_2016_Jan19/mt_keras_20190124
ARTUS_FRIENDS_TT_2016=/ceph/swozniewski/SM_Htautau/ntuples/Artus16_Nov/Artus_2016_Jan19/tt_keras_20190124
ARTUS_FRIENDS_FAKE_FACTOR_2016=/ceph/swozniewski/SM_Htautau/ntuples/Artus16_Nov/Artus_2016_Jan19/fake_factor_friends_njets_mvis_NEW_NN_Jan27
ARTUS_FRIENDS_FAKE_FACTOR_INCL_2016=$ARTUS_FRIENDS_FAKE_FACTOR_2016

# Samples Run2017
ARTUS_OUTPUTS_2017=/ceph/swozniewski/SM_Htautau/ntuples/Artus17_Jan/merged_fixedjets
ARTUS_FRIENDS_ET_2017=/ceph/swozniewski/SM_Htautau/ntuples/Artus17_Jan/et_keras_20190125
ARTUS_FRIENDS_MT_2017=/ceph/swozniewski/SM_Htautau/ntuples/Artus17_Jan/mt_keras_20190125
ARTUS_FRIENDS_TT_2017=/ceph/swozniewski/SM_Htautau/ntuples/Artus17_Jan/tt_keras_20190125
ARTUS_FRIENDS_FAKE_FACTOR_2017=/ceph/swozniewski/SM_Htautau/ntuples/Artus17_Jan/fake_factor_friends_njets_mvis_NEW_NN_Jan26
ARTUS_FRIENDS_FAKE_FACTOR_INCL_2017=$ARTUS_FRIENDS_FAKE_FACTOR_2017

# Samples Run2018
ARTUS_OUTPUTS_2018=/portal/ekpbms2/home/jbechtel/no_ee_noise_files/

# Error-handling
if [[ $ERA == *"2016"* ]]
then
    ARTUS_OUTPUTS=$ARTUS_OUTPUTS_2016
    ARTUS_FRIENDS_ET=$ARTUS_FRIENDS_ET_2016
    ARTUS_FRIENDS_MT=$ARTUS_FRIENDS_MT_2016
    ARTUS_FRIENDS_TT=$ARTUS_FRIENDS_TT_2016
    ARTUS_FRIENDS_FAKE_FACTOR=$ARTUS_FRIENDS_FAKE_FACTOR_2016
    ARTUS_FRIENDS_FAKE_FACTOR_INCL=$ARTUS_FRIENDS_FAKE_FACTOR_INCL_2016
elif [[ $ERA == *"2017"* ]]
then
    ARTUS_OUTPUTS=$ARTUS_OUTPUTS_2017
    ARTUS_FRIENDS_ET=$ARTUS_FRIENDS_ET_2017
    ARTUS_FRIENDS_MT=$ARTUS_FRIENDS_MT_2017
    ARTUS_FRIENDS_TT=$ARTUS_FRIENDS_TT_2017
    ARTUS_FRIENDS_FAKE_FACTOR=$ARTUS_FRIENDS_FAKE_FACTOR_2017
    ARTUS_FRIENDS_FAKE_FACTOR_INCL=$ARTUS_FRIENDS_FAKE_FACTOR_INCL_2017
elif [[ $ERA == *"2018"* ]]
then
    ARTUS_OUTPUTS=$ARTUS_OUTPUTS_2018
fi

# Kappa database
KAPPA_DATABASE=/portal/ekpbms1/home/jbechtel/postprocessing/2018/sm-htt-analysis/datasets.json

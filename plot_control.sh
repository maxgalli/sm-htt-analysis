source utils/setup_cvmfs_sft.sh
source utils/setup_python.sh

ERA=$1
CHANNEL=$2
VARIABLE=$3

#./shapes/produce_shapes_variables_nominal.sh $ERA $CHANNEL $VARIABLE

plotting/plot_shapes_control.py -l  --era Run${ERA} --input ${ERA}_${CHANNEL}_${VARIABLE}_shapes.root --variables $VARIABLE --channels $CHANNEL #--embedding --fake-factor

#v="NNrecoil_pt,nnmet,mt_1_nn,mt_2_nn,mt_tot_nn,pt_tt_nn,pZetaNNMissVis,pt_ttjj_nn,mTdileptonMET_nn"
#ch="mt,et,tt,em"
#plotting/plot_shapes_control.py -l  --era Run2017 --input shapes.root --variables ${v} --channels ${ch} --embedding --fake-factor
#plotting/plot_shapes_control.py -l  --era Run2017 --input shapes.root --variables ${v} --channels ${ch} --fake-factor
#plotting/plot_shapes_control.py -l  --era Run2017 --input shapes.root --variables ${v} --channels ${ch} --embedding
#plotting/plot_shapes_control.py -l  --era Run2017 --input shapes.root --variables ${v} --channels ${ch}

v="m_vis,m_vis_high,ptvis,met,puppimet,metParToZ,metPerpToZ,puppimetParToZ,puppimetPerpToZ,pt_1,pt_2,eta_1,eta_2,njets,jpt_1,jpt_2,jeta_1,jeta_2"
ch="mm"
plotting/plot_shapes_control.py -l  --era Run2017 --input shapes.root --variables ${v} --channels ${ch}
plotting/plot_shapes_control.py -l  --era Run2017 --input shapes.root --variables ${v} --channels ${ch} --category-postfix "peak"

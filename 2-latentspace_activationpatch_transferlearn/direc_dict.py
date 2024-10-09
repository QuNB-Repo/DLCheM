#Dictionary of preset directories/filepaths
direc_dict_train = {'qm9Oembsvselectrondensity':'../data/embs/embsdensityince/20240907oxygenfullspectrumwfullP/OembsVSfullocc250to5000NNmodel1/OembsVSfullocc.csv',
              'qm9Cembsvselectrondensity':'../data/embs/embsdensityince/20240910carbonwtotalED/CembsVSfullocc250to5000NNmodel1/CembsVSfullocc.csv',
              'iupacpkaembsvspka':'../data/embs/embspKA-stef/1layerNN-1L1NLbias258/601Oembs.csv',
              'qm9embeddingsvsae':'../data/embs/old/model1-10000/layer6/Cembslayer6-addedfromprevlayer.csv',
              'qm9embeddingsvsnmr': '../data/embs/embsQM9NMR/Cembslayer6-nmrgas.csv',
              'qm9embeddingswcarbacid@6000':'../data/perturbation/targetall-useavg-subavg/qm9-5000-6000-all/mlp/Caewcarbacid/carbacidd1-p5.csv'}  #X_test = X_test[:-2] NOTE missing last two carbons got deleted, to match it with 10000 QM9 molecules 




#Dictionary of preset directories/filepaths
direc_dict_test = {'qm9Oembsvselectrondensity':'../data/embs/embsdensityince/20240907oxygenfullspectrumwfullP/OembsVSfullocc250to5000NNmodel1/OembsVSfullocc.csv',
              'qm9Cembsvselectrondensity': '../data/embs/embsdensityince/20240910carbonwtotalED/CembsVSfullocc250to5000NNmodel1/CembsVSfullocc.csv',
              'qm9ED_perts':'../data/embs/embsdensityince/pertsmodel/0.csv',
             'qm9embsvselectrondensity':'../data/embs/embsdensity-ince/Cembs250-5000woccs.csv',
             'NMRstudy_cyclopentafuran':'../data/perturbation/NMRstudy_cyclopentafuran_BEST_0-2000/3.csv',
             'sym_alkene_H2ox_reactembs_perts':'../data/fgtransform/model1/sym_alkene_H2ox/rpert1.csv',
              'sym_alkene_H2ox_prodembs_perts':'../data/fgtransform/model1/sym_alkene_H2ox/ppert1.csv',
              'sym_alkene_H2ox_reactembs':'../data/fgtransform/model1/sym_alkene_H2ox/reactembsfull.csv',
              'sym_alkene_H2ox_prodembs':'../data/fgtransform/model1/sym_alkene_H2ox/prodembsfull.csv',
              'pr_alcohol_H2ox_reactembs_perts':'../data/fgtransform/model1/pr_alcohol_H2ox/r1.csv',
              'pr_alcohol_H2ox_prodembs_perts':'../data/fgtransform/model1/pr_alcohol_H2ox/p1.csv',
              'pr_alcohol_H2ox_reactembs':'../data/fgtransform/model1/pr_alcohol_H2ox/reactembsfull.csv',
              'pr_alcohol_H2ox_prodembs':'../data/fgtransform/model1/pr_alcohol_H2ox/prodembsfull.csv',
              'alcohol_H2ox_reactembs_perts':'../data/fgtransform/model1/alcohol_H2ox/r3.csv',
              'alcohol_H2ox_prodembs_perts':'../data/fgtransform/model1/alcohol_H2ox/p4.csv',
              'alcohol_H2ox_reactembs':'../data/fgtransform/alcohol_H2ox/totalrembs_Ctarget.csv',
              'alcohol_H2ox_prodembs':'../data/fgtransform/alcohol_H2ox/totalpembs_Ctarget.csv',
              'alkane_H2ox_prodembs_perts':'../data/fgtransform/model1/alkane_H2ox/p2.csv',
              'alkane_H2ox_reactembs_perts':'../data/fgtransform/model1/alkane_H2ox/r2.csv',
              'alkane_H2ox_prodembs':'../data/fgtransform/alkane_H2ox/totalpembs_Ctarget.csv',
              'alkane_H2ox_reactembs':'../data/fgtransform/alkane_H2ox/totalrembs_Ctarget.csv',
              'qm9embeddingsvsae':'../data/embs/old/model1-10000/layer6/Cembslayer6-addedfromprevlayer.csv',
              'qm9embeddingsvsnmr': '../data/embs/embsQM9NMR/Cembslayer6-nmrgas.csv',
              'qm9embeddingswcarbacid@6000':'../data/perturbation/targetall-useavg-subavg/qm9-5000-6000-all/mlp/Caewcarbacid/carbacidd1-p5.csv'}  #X_test = X_test[:-2] NOTE missing last two carbons got deleted, to match it with 10000 QM9 molecules 


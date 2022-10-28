##DeepTCR
# #clean install
# cd /n/scratch3/users/e/eha8/Rsessions/20221022_scTfh
# virtualenv deeptcr_venv
# source deeptcr_venv/bin/activate
# pip3 install DeepTCR

#startup
import matplotlib
matplotlib.use("agg") 
import os, sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from DeepTCR.DeepTCR import DeepTCR_U
from DeepTCR.DeepTCR import DeepTCR_SS
import matplotlib.colors
os.chdir("/data/deepTCR/results")


###Unsupervised CNN on repertoires (disease-based) for clustering
# Instantiate training object
DTCRUr = DeepTCR_U('./unsupervised.rep')
#Load Data from directories
DTCRUr.Get_Data(directory='../scTfh.data',Load_Prev_Data=False,aggregate_by_aa=True,
               aa_column_beta=0,v_beta_column=1,j_beta_column=2, d_beta_column=3,
               aa_column_alpha=4, v_alpha_column=5,j_alpha_column=6,count_column=8)
            
colors={"AID":matplotlib.colors.to_rgb("#F8766D"),"m564":matplotlib.colors.to_rgb("#00BFC4")}
# #Train VAE
DTCRUr.Train_VAE(Load_Prev_Data=False,var_explained=0.99)
features = DTCRUr.features
DTCRUr.Sample_Features()
DTCRUr.sample_features
print(features.shape)
DTCRUr.sample_features.to_csv("unsupervised.rep.features.csv")
DTCRUr.HeatMap_Sequences()
plt.savefig("unsupervised.rep.heatmap.seq.png")
DTCRUr.HeatMap_Samples()
plt.savefig("unsupervised.rep.heatmap.sample.ratio.png")
#Clustering
DTCRUr.Cluster(clustering_method='phenograph',write_to_sheets=True,sample=1000) #hierarchical, dbscan; use sample to downsample
DFs = DTCRUr.Cluster_DFs
df = pd.concat(DFs,keys=list(range(0,len(DFs))),names=['cluster','rowID'])
df.to_csv("unsupervised.rep.clusters.csv")
DTCRUr.UMAP_Plot(by_cluster=True, scale=5,show_legend=False)
plt.savefig("unsupervised.rep.umap.VAE.clusters.png")
DTCRUr.UMAP_Plot(Load_Prev_Data=True,by_sample=True, scale=5,show_legend=False)
plt.savefig("unsupervised.rep.umap.VAE.samples.png")
DTCRUr.UMAP_Plot(Load_Prev_Data=True,by_class=True,freq_weight=True,scale=2500,alpha=0.5)
plt.savefig("unsupervised.rep.umap.VAE.class.png")
DTCRUr.UMAP_Plot(Load_Prev_Data=True,by_class=True,freq_weight=True,scale=2500,alpha=0.5,show_legend=False)
plt.savefig("unsupervised.rep.umap2.VAE.class.png")
DTCRUr.UMAP_Plot_Samples(scale=100)
plt.savefig("unsupervised.rep.umap.samples.png")
#export UMAP
import umap
import numpy
import pickle
umap_obj = umap.UMAP()
features_umap = umap_obj.fit_transform(DTCRUr.features)
df=numpy.c_[features_umap,DTCRUr.class_id,DTCRUr.sample_id,DTCRUr.beta_sequences,DTCRUr.alpha_sequences,DTCRUr.freq,DTCRUr.counts]
numpy.savetxt("unsup.rep.umap.csv",df,delimiter=",",fmt="%s")# save train feature numpy array to csv
pickle.dump(DTCRUr.features,open("unsup.rep.obj","wb"))
#Repertoire visualization
DTCRUr.Repertoire_Dendrogram(n_jobs=40,distance_metric='KL')
DTCRUr.Repertoire_Dendrogram(lw=6,gridsize=25,Load_Prev_Data=True,
    gaussian_sigma=0.75,dendrogram_radius=0.2,repertoire_radius=0.3,color_dict=colors) #,sample_labels=True
plt.savefig("unsupervised.rep.repertoire.dendrogram.png")
#KNN classification (testing ML predictive capability)
DTCRUr.KNN_Repertoire_Classifier(metrics=['AUC'],distance_metric='KL',plot_metrics=True,by_class=True)
# DTCRUr.KNN_Repertoire_Classifier(metrics=['AUC','Recall'], #,"Precision","F1_Score"
#   distance_metric='KL',plot_metrics=True,by_class=True,Load_Prev_Data=True) #metrics: KL, correlation, euclidean, wasserstein, JS; Sample to speed up
plt.savefig("unsupervised.rep.KNN.scores.png")
DTCRUr.KNN_Repertoire_DF.to_csv("unsupervised.rep.KNN.csv")
#Motif Identification (for WebLogo)
DTCRUr.Motif_Identification(group='m564',by_samples=True)
DTCRUr.Motif_Identification(group='AID',by_samples=True)
#Structural Diversity
DTCRUr.Structural_Diversity(sample=1000) #sample=500
print(DTCRUr.Structural_Diversity_DF)
DTCRUr.Structural_Diversity_DF.to_csv("unsupervised.rep.structural.diversity.csv")


# 
# ###Supervised repertoire classification (TILs --> predict treatment category)
# from DeepTCR.DeepTCR import DeepTCR_WF
# # Instantiate training object
# DTCR_WF = DeepTCR_WF('./supervised.rep.class')
# #Load Data from directories
# # DTCR_WF.Get_Data(directory='../tutorial.data/Rudqvist',Load_Prev_Data=False,aggregate_by_aa=True,
# #                aa_column_beta=1,count_column=2,v_beta_column=7,d_beta_column=14,j_beta_column=21)
# # DTCR_WF.Get_Data(directory='../training.data/disease2',Load_Prev_Data=False,aggregate_by_aa=True,
# #                aa_column_beta=0,v_beta_column=1,j_beta_column=2)
# DTCR_WF.Get_Data(directory='../scTfh.data',Load_Prev_Data=False,aggregate_by_aa=True,
#                aa_column_beta=0,v_beta_column=1,j_beta_column=2, d_beta_column=3,
#                aa_column_alpha=4, v_alpha_column=5,j_alpha_column=6,count_column=8)
# #train model (method 1)
# DTCR_WF.Get_Train_Valid_Test(test_size=0.25)
# DTCR_WF.Train()
# DTCR_WF.AUC_Curve()
# plt.savefig("supervised.rep.method1.AUC.samples.png")
# #train model (method 2)
# DTCR_WF.K_Fold_CrossVal(combine_train_valid=True, hinge_loss_t = 0.1,train_loss_min = 0.1,
#   folds=5,num_concepts=64, size_of_net="small")
# DTCR_WF.AUC_Curve()
# plt.savefig("supervised.rep.method2.AUC.samples.png")
# #train model (method 3 - Monte Carlo)
# DTCR_WF.Monte_Carlo_CrossVal(folds=25,LOO=4,epochs_min=10,num_concepts=64,size_of_net="small",
#                              train_loss_min=0.1,hinge_loss_t=0.1,combine_train_valid=True)
# DTCR_WF.AUC_Curve()
# DTCR_WF.AUC_Curve(figsize=(5,6))
# plt.savefig("supervised.rep.method3.AUC.samples.png")
# #representative sequences and motifs
# DTCR_WF.Representative_Sequences()
# print(DTCR_WF.Rep_Seq['m564'])
# df = pd.concat(DTCR_WF.Rep_Seq,keys=list(DTCR_WF.Rep_Seq.keys()),names=['Disease','rowID'])
# df.to_csv("supervised.rep.rep_seq.csv")
# DTCR_WF.Motif_Identification('m564')
# DTCR_WF.Motif_Identification('AID')
# # The motifs can then be found in fasta files in the results folder underneath (label)(alpha/beta)Motifs. 
# # These fasta fiels can then be used with "https://weblogo.berkeley.edu/logo.cgi" for motif visualization.
# #visualize learned latent space
# DTCR_WF.UMAP_Plot(by_class=True,freq_weight=True,scale=1000)
# DTCR_WF.UMAP_Plot(by_class=True,freq_weight=True,scale=2000,set='test') #specify if using train, valid, or test data
# DTCR_WF.UMAP_Plot(Load_Prev_Data=True,by_class=True,freq_weight=True,scale=2000,alpha=0.5)
# plt.savefig("supervised.rep.umap.png")
# #repertoire relatedness
# DTCR_WF.Repertoire_Dendrogram(gridsize=50,gaussian_sigma=0.75,lw=6,dendrogram_radius=0.3)
# DTCR_WF.Repertoire_Dendrogram(lw=6,gridsize=25,Load_Prev_Data=True,
#     gaussian_sigma=0.75,dendrogram_radius=0.2,repertoire_radius=0.3,color_dict=colors) #,sample_labels=True
# plt.savefig("supervised.rep.dendrogram.png")
# #export UMAP
# import umap
# import numpy
# import pickle
# umap_obj = umap.UMAP()
# features_umap = umap_obj.fit_transform(DTCR_WF.features)
# df=numpy.c_[features_umap,DTCR_WF.class_id,DTCR_WF.sample_id,DTCR_WF.beta_sequences,DTCR_WF.alpha_sequences,DTCR_WF.freq,DTCR_WF.counts]
# numpy.savetxt("sup.rep.umap.csv",df,delimiter=",",fmt="%s")# save train feature numpy array to csv
# pickle.dump(DTCR_WF.features,open("unsup.rep.obj","wb"))
# 
# 
# ###Supervised sequence regression (predicting a continuious rather than descrete variable from TCRs)
# import pandas as pd
# from DeepTCR.DeepTCR import DeepTCR_SS
# import numpy as np
# import os
# import matplotlib.pyplot as plt
# df = pd.read_csv('../../BMchim.counts_regression.csv')
# # df = pd.read_csv('../tutorial.data/10x_Data/Data_Regression.csv')
# #define count variable
# Y = np.log2(np.array(df['m564'])+1) 
# # Y = np.log2(np.array(df['A0201_GILGFVFTL_Flu-MP_Influenza'])+1) 
# #Instantiate object
# DTCR = DeepTCR_SS('./supervised.regression')
# DTCR.Load_Data(alpha_sequences=np.array(df['alpha']),beta_sequences=np.array(df['beta']),Y=Y)
# #train (method #1)
# DTCR.K_Fold_CrossVal(folds=5) 
# #train (method #2)
# DTCR.Monte_Carlo_CrossVal(folds=5)
# #density distirbution of predicted vs ground truth
# import matplotlib.pyplot as plt
# from scipy.stats import gaussian_kde
# def Plot_Func(x,y):
#     xy = np.vstack([x, y])
#     z = gaussian_kde(xy)(xy)
#     r = np.argsort(z)
#     x ,y, z = x[r], y[r], z[r]
#     plt.figure()
#     plt.scatter(x,y,s=15,c=z,cmap=plt.cm.jet)
#     plt.xlim([0,10])
#     plt.ylim([0,10])
#     plt.xlabel('Predicted',fontsize=12)
#     plt.ylabel('Log2(counts)',fontsize=12)
#     
# Plot_Func(np.squeeze(DTCR.y_pred),np.squeeze(DTCR.y_test))
# plt.savefig("supervised.seq.reg.png")
# Plot_Func(np.squeeze(DTCR.predicted),Y) #if Monte Carlo
# plt.savefig("supervised.seq.reg2.png")
# #Representative sequences
# DTCR.Representative_Sequences(top_seq=25)
# DTCR.Rep_Seq['High']
# df = pd.concat(DTCR.Rep_Seq,keys=list(DTCR.Rep_Seq.keys()),names=['Level','rowID'])
# df.to_csv("supervised.rep.reg.rep_seq.csv")
# 
# 
#  




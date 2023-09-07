### this is the first part of making heatmaps using features selected by the significant latent factors using ER. This code orders the features ####

## this code is the first part of making heatmaps for Schisto univariate features
import pandas as pd
import numpy as np

Xdf=pd.read_csv("/ix/djishnu/Trirupa/Schisto_Proj/ERinputs/MFI50/Xc2_c3_AbScAg_MFI50_30May23.csv", sep=",")
Ydf=pd.read_csv("/ix/djishnu/Trirupa/Schisto_Proj/ERinputs/Yc2_c3_AbScAg_30May23.csv", sep=",")
Xdf=Xdf.drop(columns=["Unnamed: 0"])
Ydf=Ydf.drop(columns=["Unnamed: 0"])
XY_df=pd.concat([Xdf,Ydf],axis=1)

## SEA+Egg- vs SEA+Egg+
Z1=["IgG1.SEA", "IgG1.Gst.SM25", "IgG3.SEA", "IgG3.Calumenin", "IgA.Gst.SM25", "IgA1.SEA", "IgA1.MEG", "IgA2.MEG", "IgA2.Calumenin", "FcR2A.Calumenin","FcR3B.Calumenin", "SNA.MEG", "SNA.CD63"]
Z2=["IgG2.Gst.SM25", "IgG2.Calumenin", "IgG3.CD63", "IgG4.SEA", "IgG4.CD63", "FcREpsilonI.SEA", "FcR2A.Calumenin"]
Z3=["IgG.Calumenin", "IgA.Calumenin", "IgA1.SEA", "IgA1.Gst.SM25", "IgA1.MEG", "IgA1.CD63", "IgA1.Calumenin", "IgA2.Calumenin"]
Z4=["IgG.Calumenin", "IgG1.Calumenin", "FcR1.CD63", "FcR1.Calumenin", "FcR3B.Calumenin"]

def Union(lst1,lst2,lst3,lst4): ### adjust according to number of Zs
    combined_list = lst1 + lst2 +lst3+lst4
    uniq_set= set(combined_list)
    final_list=list(uniq_set)
    return final_list

all_feat=Union(Z1,Z2,Z3,Z4)

XY_df = XY_df.rename(columns = {'UnsupY':'Y'})
Y1_df=XY_df.loc[XY_df["Y"]==1]
Y0_df=XY_df.loc[XY_df["Y"]==0]


mean_diff={}
for feat in all_feat:
    mean_0=Y0_df[feat].mean(axis=0)
    mean_1=Y1_df[feat].mean(axis=0)
    mean_diff[feat]= [mean_0,mean_1]
mean_df=pd.DataFrame.from_dict(mean_diff, orient='index', columns=['mean_0', 'mean_1'])
mean_df.reset_index(inplace=True)

## calculation mean(log2(fold change)) ##
mean_df = mean_df.rename(columns = {'index':'Feature'})
mean_df["FC"]=abs((mean_df['mean_1']-mean_df["mean_0"])/mean_df["mean_0"])
mean_df["log2FC"]=np.log2(mean_df["FC"])
meandf_sorted = mean_df.sort_values(by='log2FC', ignore_index=True,na_position='first')
feat_sorted= list(meandf_sorted['Feature'])

XYsort_df =XY_df.reindex(columns=feat_sorted) 
XYsort_df.to_csv("/ix/djishnu/Trirupa/Schisto_Proj/forHeatmap/data/c2_c3_featsorted_spec0.1_patid_X.csv", sep=",",index=None)
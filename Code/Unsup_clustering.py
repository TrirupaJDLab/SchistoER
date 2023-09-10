### this code performs unsupervised clustering on the SEA+Egg- group of samples ####

## loading libraries 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
%matplotlib inline 
from sklearn.cluster import KMeans
from sklearn import datasets
from scipy.spatial.distance import cdist
from yellowbrick.cluster import KElbowVisualizer
from sklearn.cluster import KMeans


dataset=pd.read_csv("Class2_XY.csv", sep=",")
AB_data=dataset.drop(columns=["Y"])

## identifying the number of clusters in the dataset
model = KMeans()
# k is range of number of clusters.
visualizer = KElbowVisualizer(model, k=(2,30),metric='silhouette', timings= True)
visualizer.fit(AB_data)        # Fit the data to the visualizer
visualizer.show(outpath="kmeans_silhouetteIndex.pdf") 

## generating outcome labels based on kmeans clustering
kmeans = KMeans(n_clusters=2, random_state=0) 
labels=kmeans.fit(AB_data)
AB_data["labels"]=labels.labels_
AB_data.to_csv("class2_withunsupLabel_patId.csv",sep=",", index=None)

import cooler
import numpy as np
import pandas as pd
import glob, os
import os.path
import matplotlib.pyplot as plt
from sklearn import decomposition
import matplotlib


######
name_list=['FA-DpnII','FA+DSG_DpnII','FA+DSG-MNase']

cmap =['r','b','g']

dot_list=list1 # (# of dots for FA-DpnII, # of dots for FA+DSG-DpnII, # of dots for FA+DSG-MNase)
anchor_list=list2 # (# of dots_anchors for FA-DpnII, # of dots_anchors for FA+DSG-DpnII, # of dots_anchors for FA+DSG-MNase)

fig, ax = plt.subplots(ncols=1,figsize=(8,8))
for i in range(10):
	ax.scatter(dot_list[i],anchor_list[i],color=cmap1[i], s=15, marker="o")
	ax.annotate(name_list[i], (dot_list[i],anchor_list[i]),fontsize=9,ha='right',color=cmap[i])
lo, hi = 0, np.max(anchor_list)
#plt.plot([2*lo, 2*hi], [lo, hi],'k-')
plt.plot([lo, hi], [2*lo, 2*hi],'r')
ax.set_xlabel("# of dots",fontsize=15)
ax.set_ylabel("# of dot anchors",fontsize=15)
plt.savefig("dots_vs_anchors.pdf")
#plt.show()
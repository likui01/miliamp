# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 22:41:51 2016

@author: likui
"""

import csv
import numpy as np
import sys
import matplotlib.pyplot as plt #importing plotting module
from Bio.Data import IUPACData as ip
from sklearn import svm
from sklearn.metrics import roc_curve, auc

from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
from roc import *
from sklearn import cross_validation
from sklearn.cross_validation import KFold
from sklearn.grid_search import GridSearchCV
from itertools import product
from sklearn.cross_validation import StratifiedKFold as SK
from roc import *
from numpy.random import choice
import Bio
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import numpy as np
import matplotlib.pyplot as plt #importing plotting module
from Bio import SeqIO
from temp2 import kmer, encode1,encode2
handle = open("S3.fasta", "rU")
records = list(SeqIO.parse(handle, "fasta"))
handle.close()
seq=[]
target=[]
for i in range(34):
    seq.append(str(records[i].seq))
    target.append(np.zeros((1,len(seq[i]))))
    

target[0][0,16:32]=1
target[0][0,87:99]=1    
target[1][0,12:28]=1    
target[2][0,7:20]=1
target[2][0,13:21]=1
target[2][0,19:29]=1
target[2][0,20:30]=1
target[2][0,29:37]=1
target[3][0,0:93]=1
target[4][0,59:70]=1
target[5][0,10:25]=1
target[5][0,24:35]=1
target[5][0,29:40]=1
target[5][0,36:42]=1
target[6][0,20:31]=1
target[6][0,32:41]=1
target[6][0,58:71]=1
target[6][0,82:89]=1
target[6][0,90:96]=1
target[7][0,10:20]=1
target[7][0,100:110]=1
target[7][0,115:126]=1
target[7][0,145:152]=1
target[8][0,14:20]=1
target[9][0,80:125]=1
target[10][0,0:35]=1
target[10][0,35:67]=1
target[11][0,97:103]=1
target[12][0,0:29]=1
target[12][0,100:118]=1
target[13][0,172:230]=1
target[14][0,217:289]=1
target[15][0,12:18]=1
target[16][0,10:17]=1
target[17][0,491:509]=1
target[18][0,537:545]=1
target[19][0,8:34]=1
target[20][0,4:14]=1
target[20][0,24:34]=1
target[20][0,55:61]=1
target[21][0,83:105]=1
target[21][0,104:125]=1
target[21][0,147:153]=1
target[21][0,153:163]=1
target[21][0,155:171]=1
target[21][0,179:196]=1
target[21][0,208:231]=1
target[22][0,31:41]=1
target[22][0,41:50]=1
target[23][0,65:72]=1
target[24][0,111:157]=1
target[25][0,6:21]=1
target[25][0,19:34]=1
target[25][0,42:57]=1
target[26][0,23:32]=1
target[27][0,0:142]=1
target[28][0,0:0:12]=1
target[29][0,5:12]=1
target[30][0,34:44]=1
target[30][0,48:59]=1
target[30][0,59:68]=1
target[30][0,68:82]=1
target[30][0,85:95]=1
target[31][0,588:600]=1
target[32][0,9:20]=1
target[32][0,104:115]=1
target[33][0,0:88]=1



    
        
protein=[]
labels=[]
for i in range(34):
    f=open('pro_all'+str(i)+'.txt','w')
    ww=[]#ds=np.array(3)# to store the score + labels 
    w1=[]
    for yy in range(len(seq[i])-6):
        aa=seq[i][yy:yy+6]
        bb=sum(target[i][0,yy:yy+6])
       # ww=ww+[aa]
        if bb>1.0:
            w1.append(1)
            f.write(aa)
            f.write(","+str(1))  
            f.write("\n")  
        elif bb==0.0:
            f.write(aa)
            w1.append(0)
            f.write(","+str(0)) 
            f.write("\n")  
        else:
            pass
    f.close()     
       # f.write("\n")    
      #  protein.append(ww)   
       # labels.append(w1)


#test=protein[0:12]
#test_labels=labels[0:12]
#train=protein[12::]
#train_labels=labels[12::]
#score=[]
#tlb=[] 
#for i in range(1):
#    test=protein[i]
#    test_labels=labels[i]
#    train=protein[0:i]+protein[i+1::]
#    train_labels=labels[0:i]+labels[i+1::]
#    Xtr=[]
#    Xla=[]
#    for ii in range(33):
#        for yy in range(len(train[ii])):
#           Xtr.append(train[ii][yy])
#           Xla.append(train_labels[ii][yy])
#           
#    scr=svm.SVC(kernel='linear', class_weight='auto')
#    parameters={'C':[0.00001,0.00001,0.0001,0.001,0.01,0.1,1,10]}#,'gamma':[0.0001,0.001,0.01,0.1,1,10,100,1000]}
#    clf = GridSearchCV(scr, parameters, cv=2, scoring='roc_auc') 
#    clf.fit(Xtr,Xla)
#       # ad=clf.best_params_
#    print clf.best_params_
#    print clf.best_score_         
#    clff= clf.best_estimator_# initaiate the SVM instance 
#    clff.fit(Xtr, Xla)        
#           
#    
#    
#    for jj in range(len(test)):      
#        qq=clff.decision_function(test[jj])
#        score.append(qq)
#        tlb.append(test_labels[jj])
#
#fpr, tpr, thresholds = roc_curve(tlb,score) # plotting ROC 
#a=auc(fpr,tpr) 
#print a           
#plt.plot(fpr,tpr,marker='.')
#plt.xlabel('fpr')
#plt.ylabel('tpr')
#plt.xlim([0,1])
#plt.ylim([0,2])
# # plt.title('ROC Linear kernel for SVM ')
#plt.grid()
#
#
#plt.legend(['auc= '+str(round(a,3))],loc=4)





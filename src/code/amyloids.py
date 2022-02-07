# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 00:57:04 2016

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
from roc import*
from features import encode2, kmer,weightmer, BOL 

def read_33pro():
    amino=[]
    target=[]
    ai=[]
    lab=[]
    file=open('33skiplist.txt',"r")
    lines=file.readlines()
    file.close()
    for i in range(len(lines)):
        if (lines[i][0] == '>' and i>0) :
            amino.append(ai)
            target.append(lab)
            ai=[]
            lab=[]
        else:  
            if(lines[i][0] != '>'):
                n=lines[i].replace('\t',',').split(',')
                ai.append(n[0])
                lab.append(int(n[1][0]))        
    return amino, target        
        
    

def read_data():
    n=[]
    amino=[]
    target=[]
   
    amino1=[]
    target1=[]
    n1=[]
     ###reading training data 
    file = open('S1.txt',"r") #open data file
    lines = file.readlines() # readlines
    file.close()
    for i in range(len(lines)):
             n.append(lines[i].replace('\t',',').split(','))  # remove tabs 
             amino.append(n[i][0]) # save sequence in a list 
             target.append(int(n[i][1])) # save labels in list
    ####reading testing data  
           
    file = open('S2.txt',"r") 
    lines1 = file.readlines() 
    file.close()
    for ii in range(len(lines1)-1):
             n1.append(lines1[ii].replace('\t',',').split(','))  
             amino1.append(n1[ii][0]) 
             target1.append(int(n1[ii][1]))     
    return amino, target, amino1, target1         
def Compute_features (amino):  
    Xtr=[]      
    for i in range(len(amino)):
        Xtr.append(kmer(amino[i])) # for the taining sequences count 1 mer features 
    return Xtr    
# details about SVM function  http://scikit-learn.org/stable/modules/generated/sklearn.svm.SVC.html    
if __name__ == '__main__':  
    amino, target, amino1, target1=read_data()
    Xtr=Compute_features(amino)
    scr=svm.SVC(C=0.1,kernel='linear', class_weight='auto')         
    scr.fit(Xtr, target)
    dd=[] # calculated targets 
    score=[]
    prob=[]
    for kk in range(len(amino1)): # use sliding window to calculate  the score 
        ww=[]#ds=np.array(3)# to store the score + labels 
        for yy in range(len(amino1[kk])-6):
            aa=kmer(amino1[kk][yy:yy+6])
            ww.append(np.hstack([scr.decision_function(aa)[0][0],scr.predict(aa)])) # using SVM to predict  the label 
        ww=np.array(ww)
        ind=np.argmax(ww[:,0])  
        #seq=kmer(amino1[kk][ind:ind+6])
        #prob.append(clf.predict_proba(seq)[0][1])
        prob.append(ww[ind,0])
        dd.append(ww[ind,1])

  
         
         
    fpr, tpr, thresholds = roc_curve(target1,prob) # plotting ROC 
    a=auc(fpr,tpr)
    print a
    precision, recall, _ = precision_recall_curve(target1,prob)
    a= average_precision_score(target1, prob)
    print a
    #plt.plot(precision, recall, marker='o')
    ###
    #plt.xlabel('recall')
    #plt.ylabel('precision')
    ##plt.title('PR curve  ')
    line1=plt.plot(fpr,tpr, marker=',')
    #plt.title('PR rbf kernel of SVM')
    plt.xlabel('fpr')
    plt.ylabel('tpr')         
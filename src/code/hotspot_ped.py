# -*- coding: utf-8 -*-
"""
Created on Thu May 31 01:19:57 2018

@author: farzeen
"""# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 01:36:58 2016

@author: likui
"""

from CDhitMIL import *
from amyloids import *
from customized import *
import numpy as np
import pickle 
from sklearn.externals import joblib
def linear_svm(): # testing on DS2, training on DS1
    amino, target, amino1, target1=read_data()
    amino33,target33=read_33pro()
    Xtr=Compute_features(amino)
    scr=svm.SVC(C=.1,kernel='linear', class_weight='auto')         
    scr.fit(Xtr, target)
    dd=[] # calculated targets 
    
    let=[]
    see=[]
    for j in range(33):
        prob=[]
        for kk in range(len(amino33[j])): 
            aa=kmer(amino33[j][kk])
            prob.append(scr.decision_function(aa)[0][0])
            dd.append(scr.predict(aa)[0])
        max_score=np.max(prob)
        for jj in range(len(prob)):
            if(target33[j][jj]==1):
                prob[jj]=max_score
        see=see+prob
        let=let+target33[j]
    return see,let 

def Linear_MIL(): # testing on DS2 and training on DS1+DS2
    bags_ds1=create_bags('S1.txt')
    bags_ds2=create_bags('S2.txt')
   
    bags=bags_ds1+bags_ds2

    classifier=linclass(epochs=5000, Lambda=0.0001)

    Folds=create_folds(bags,5)
    prob=[]
    target=[]
    for i in range(5):
        train=[bags[Folds[i].train_bags[ii]] for ii in range(len(Folds[i].train_bags))]
        test=[bags[Folds[i].test_bags[ii]] for ii in range(len(Folds[i].test_bags))]
        classifier.train(train,merge=True)
        for i in range(34):
             bag3=create_bags('pro_all'+str(i)+'.txt')
             score=classifier.test(bag3)
             labels=[]
###        
             labels+=[b.label for b in bag3]  
             max_score=np.max(score)
             for jj in range(len(score)):
                 if(labels[jj]==1):
                     score[jj]=max_score
             prob.append(score)
             target.append(labels)
        

    sco=[]
    lab=[]    
    for i in range(5):
      for k in range(len(prob[i])):
          sco.append(prob[i][k])
          lab.append(target[i][k])
    return lab,sco
def Customized_classifier(): #testing on DS2, training on DS1+DS2+Muatants-aggrescan 
       test1=create_bags('S1.txt')
       test2=create_bags('S2.txt')
      
#       train1,folds=create_mutants()
       bag_aggre, bag_other=sep_aggre()
       bags=test1+test2

       classifier=linclassrank(epochs=5000,Lambda=0.001)#0.0001
       Folds=create_folds(bags,5)
       
       prob=[]
       target=[]
       
       for i in range(5):
                train=[bags[Folds[i].train_bags[ii]] for ii in range(len(Folds[i].train_bags))]
                test=[bags[Folds[i].test_bags[ii]] for ii in range(len(Folds[i].test_bags))]
#                classifier=linclassrank(epochs=5000,Lambda=1000)#0.0001
                classifier.train(train,bag_other)
                for i in range(34):
                 bag3=create_bags('pro_all'+str(i)+'.txt')
                 score=classifier.test(bag3)
                 labels=[]
    ###        
                 labels+=[b.label for b in bag3]  
                 max_score=np.max(score)
                 for jj in range(len(score)):
                     if(labels[jj]==1):
                         score[jj]=max_score
                 prob.append(score)
                 target.append(labels)
       joblib.dump(classifier,'hot_spot_MIL_rank.pkl')
       sco=[]
       lab=[]    
       for i in range(5):
          for k in range(len(prob[i])):
              sco.append(prob[i][k])
              lab.append(target[i][k])    
       return lab,sco       
def readscore(filename):
     sco=[]
     file = open(filename,"r") #open data file
     lines = file.readlines() # readlines
     file.close()    
     n=[]
     for i in range(len(lines)):
             n.append(lines[i].replace('\n',',').split(','))  # remove tabs 
             sco.append(float(n[i][0])) # save score 
             
     return sco
def readwaltz():
     sco=[]
     lab=[]
     file = open('waltzscore.txt',"r") #open data file
     lines = file.readlines() # readlines
     file.close()    
     n=[]
     for i in range(len(lines)):
             n.append(lines[i].replace('\t',',').split(','))  # remove tabs 
             sco.append(float(n[i][0])) # save score 
             lab.append(float(n[i][1]))
     return lab,sco  
if __name__ == '__main__':
#    lab1,prob=linear_svm() # simple SVM 
#    fpr, tpr, thresholds = roc_curve(prob,lab1)
#    plt.plot(fpr,tpr)
#    a=auc(fpr,tpr)
#    print(a)
#    lab2,sco1=Linear_MIL() # MIL
#    fpr1, tpr1, thresholds = roc_curve(lab2,sco1)
#    plt.plot(fpr1,tpr1)
#    a1=auc(fpr1,tpr1)
#    print(a1)
#    precision, recall, _ = precision_recall_curve(lab2,sco1)
#    ass= average_precision_score(lab2, sco1)
#    print ass
#    plt.plot(precision, recall, marker='o')
#    precision1, recall1, _ = precision_recall_curve(prob,lab1)
#    ass2= average_precision_score(prob, lab1)
#    print ass2
#    plt.plot(precision1, recall1, marker='o')
    lab3,sco2=Customized_classifier() #linclassrank classifier
    fpr, tpr, thresholds = roc_curve(lab3,sco2)
    plt.plot(fpr,tpr)
    a=auc(fpr,tpr)
    print(a)
#    lab4,sco3=readwaltz()
#    sco4=readscore('mettamyscore.txt')
#    sco5=readscore('agreescanmax.txt')
#    sco6=readscore('appnnscore.txt')
#    aa=[]
#    for i in [(lab1,prob,'.','-'),(lab2,sco1,'o','-'),(lab3,sco2,"3",'-.'),(lab1,sco4,'>','-'),(lab1,sco5,'v',':'),(lab1,sco6,'2','--'),(lab4,sco3,'+',':')]:
#        
#        fpr, tpr, thresholds = roc_curve(i[0],i[1]) # plotting ROC 
#        plt.plot(fpr,tpr, marker=i[2],linestyle=i[3])
#        a=auc(fpr,tpr)
#        aa.append(a)
#
#    plt.xlabel('Fpr')
#    plt.ylabel('Tpr')         
#    plt.xlim([-0.1,0.5])
#
#    plt.legend(['Linear SVM:'+str(round(aa[0],3)*100),'MIL:'+str(round(aa[1],3)*100),'MIL-Rank:'+str(round(aa[2],3)*100),'MetAmyl:'+str(round(aa[3],3)*100),'Aggrescan:'+str(round(aa[4],3)*100),'APPNN:'+str(round(aa[5],3)*100),'Waltz:'+str(round(aa[6],3)*100)],loc=4)
#    plt.grid()








# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 00:00:55 2016

@author: likui
"""

from customized import *
from CDhitMIL import *
from sklearn.externals import joblib

def MIL_Rank():
   
    test1=create_bags('S1.txt')
    test2=create_bags('S2.txt')
    train1,folds=create_mutants()
    bag_aggre, bag_other=sep_aggre()
    
    bags=test1+test2
    classifier=linclassrank(epochs=5000,Lambda=0.0001)
       
    classifier.train(bags,bag_other,Beta=0.5)
    score=classifier.test(bag_aggre)
    joblib.dump(classifier,'mute_pred.pkl')
    labels=[]
    ids=[]
    #    ###        
    labels+=[b.label for b in bag_aggre]    
    ids+=[b.id for b in bag_aggre]
   
    return score, labels
def aggrescab():
    aggre=[-16,15,29,5,-68,-63,16,-118,-15,-15,-62,-49,-12,-10,32,2,2,2,1,2,2,-1,-5,-3,9,17,11,21,-59,16,-61,-23,-106]
    tar=[-1,1,1,1,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,-1,-1,1,1,1,1,-1,1,-1,-1,1]        
    return aggre,tar
def MIL():
   
    test1=create_bags('S1.txt')
    test2=create_bags('S2.txt')
    train1,folds=create_mutants()
    bag_aggre, bag_other=sep_aggre()
    bags=test1+test2
    classifier=linclassrank(epochs=5000,Lambda=0.00001)
       
    classifier.train(bags,bag_other,Beta=0.0)
    score=classifier.test(bag_aggre)
    labels=[]
    ids=[]
    #    ###        
    labels+=[b.label for b in bag_aggre]    
    ids+=[b.id for b in bag_aggre]

    return score,labels
    
if __name__ == '__main__':
   sco,lab=MIL_Rank()
   sco1,lab1=aggrescab()
   sco2,lab2=MIL()
   aa=[]
   for i in [(lab,sco,'.','-'),(lab1,sco1,'o','-'),(lab2,sco2,"3",'-.')]:
        fpr, tpr, thresholds = roc_curve(i[0],i[1]) # plotting ROC 
        plt.plot(fpr,tpr, marker=i[2],linestyle=i[3])
        a=auc(fpr,tpr)
        aa.append(a)

   plt.xlabel('Fpr')
   plt.ylabel('Tpr')         
   plt.xlim([-0.1,1])

   plt.legend(['MIL_Rank:'+str(round(aa[0],3)*100),'Aggrescan:'+str(round(aa[1],3)*100),'MIL:'+str(round(aa[2],3)*100)],loc=4)
   plt.grid()
   plt.show() 
    
    
    
    
    
    
    
    
    
    
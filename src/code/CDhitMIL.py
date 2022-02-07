# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 22:20:12 2
to test sequences train on DS1+DS2 by 5 fold cv created by CD-hit at 40% similarity
MIL linear and Non-linear
@author: likui
"""

from classifiers import *
import numpy as np
from scipy.io import *
from bag import *
from llc import *
from cv import *
from classifiers import *
from roc import*
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
import matplotlib.pyplot as plt #importing plotting module
#from roc import *
from temp2 import *
from sklearn.cross_validation import StratifiedKFold as SK
def create_folds(bags, no_of_folds):
    """
    Creates folds from the given data.
    Takes a list of bags and the desired number of folds as input.
    Returns a list of fold objects
    """
    file = open('db2novel.txt',"r") 
    lines1 = file.readlines() 
    file.close()
    i=0
    lin=[]
    while (i<(len(lines1))):
        
        clus=[]
        if (lines1[i][0]=='>'):
            i=i+1
            
            while (lines1[i][0]!='>'):
                n=(lines1[i].replace('\t',',').split(','))
                n11=(n[2].replace('\t',',').replace('|',',').split(','))
                clus.append(int(n11[1]))
                i=i+1
                if (i==len(lines1)):
                    break
       
        lin=lin+[clus]
    
   
    ac=[]
    aa=lin[0]+lin[2]
    ab=lin[1]+lin[3]+lin[4]
    for i in range(5,15):
        ac=ac+lin[i]
    ad=[]    
    for i in range(15,40): 
        ad=ad+lin[i]
    ae=[]    
    for i in range(40,92):
        ae=ae+lin[i]
    testfolds=[aa]+[ab]+[ac]+[ad]+[ae]   
    foldstrain=[ab+ac+ad+ae]+[aa+ac+ad+ae]+[aa+ab+ad+ae]+[aa+ab+ac+ae]+[aa+ab+ac+ad]
    folds=[fold() for _ in range(no_of_folds)]
#    target1=[bags[i].label for i in range(304,len(bags))]
#    Kf=SK(target1, n_folds=no_of_folds, shuffle=False,random_state=None) 
    train_index=np.array(foldstrain)
    test_index=np.array(testfolds)
     
    for ii in range(5):
            folds[ii].train_bags=range(304)
            folds[ii].train_bags=np.append(folds[ii].train_bags,np.array(train_index[ii])+304)
            folds[ii].test_bags=np.array(test_index[ii])+304
            

    return folds

def compute_gammas(bags, C=50, beta=0.1, k=50):
    data_points=None
#    for bag_i in (bags):
#        
#        if issparse(bags[0].instances):
#            for example_j in (bag_i.instances).toarray():
#                data_points+=[example_j]
#        else:
        #for example_j in (bag_i.instances):
#    import pdb; pdb.set_trace() 
    for b in bags:
        b.gammas=[]
        for ins in b.instances:
            if data_points==None:
                data_points=np.array(ins[np.newaxis,:])
            else:
                data_points=np.vstack((data_points, np.array(ins[np.newaxis,:])))
            #data_points+=[example_j] 
   # data_points=vstack(data_points)
#    import pdb; pdb.set_trace()
    G, C, Y=llc(X=data_points,C=C, beta=beta, k=k)
    
    gamma_index=0
    for bag_index in range(len(bags)):
        for ins_index in range(np.shape(bags[bag_index].instances)[0]):
            bags[bag_index].gammas+=[G[gamma_index]]
            gamma_index+=1
        bags[bag_index].gammas=np.array(bags[bag_index].gammas)
    if issparse(bags[0].instances):        
        for bag in bags:
            bag.gammas=lil_matrix(bag.gammas).tocsr()

    return G, C, Y  
def create_bags(filename):
    n=[]    
    seq=[]
    target=[]
    file = open(filename,"r") 
    lines = file.readlines() 
    file.close()    
    for i in range(len(lines)):
         n.append(lines[i].replace('\t',',').split(','))  # remove tabs 
         seq.append(n[i][0]) # save sequence in a list 
         target.append(int(n[i][1])) # save labels in list
         
    target = 2 * np.array(target) - 1  
    if len(seq[0])==6:
         feat1=[]
         for i in range(len(seq)):
             aa=np.empty((1,20))
             ab=np.hstack((kmer(seq[i])))
             aa=np.matrix(ab)
             #print np.shape(aa.T)
             feat1.append(aa)
      
         bags=[Bag() for i in range(len(feat1))]
         for b in range(len(bags)):
            
            bags[b].instances=np.array(feat1[b])
            bags[b].label=target[b]
            bags[b].peta=[1.0]
             
    else:
#        import pdb; pdb.set_trace()     
        feat1=[]
        for kk1 in range(len(seq)): 
            ww1=[]
            for yy in range(len(seq[kk1])-6):
                aa=np.empty((1,20))
                aa=np.array(np.hstack((kmer(seq[kk1][yy:yy+6]))))
               # print type(aa)
                ww1=ww1+[aa]
#            print np.shape(ww1) 
            feat1.append(np.array(ww1))  
    
        
        bags=[Bag() for i in range(len(feat1))]
        for b in range(len(bags)):
           # import pdb; pdb.set_trace()
            for ii in range(len(feat1[b])):
                bags[b].addInstance(feat1[b][ii])
            bags[b].label=target[b]
            bags[b].peta=[1.0]*(np.shape(feat1[b])[0])        
    return bags    
            
def readBags(filename_f, filename_l, **kwargs):
    feature_matrix=loadmat(filename_f)['arrr'][0]
    labels=loadmat(filename_l)['label'][0]
    bags=[]
    for i in range(len(feature_matrix)):
        f=np.array(feature_matrix[i].T)
       # f=np.array(np.hstack((feature_matrix[i][:2].T, feature_matrix[i][2:].T)))
        #import pdb; pdb.set_trace()
        b=Bag()
        b.instances=f
        b.label=labels[i]
        b.peta=[1.0]*(np.shape(f)[0])
        bags.append(b)
    return bags
def readtestinstance(filename):
    protein=[]
    instances=[]
    target=[]
    labels=[]

     ###reading instances  
    file = open(filename,"r") #open data file
    lines = file.readlines() # readlines
    file.close()
    i=0
    while i < len(lines):
          if lines[i][0]=='>':
              if not len(instances):
                  pass
              else:
                  protein.append(instances )
                  labels.append(target)
                 # print len(instances)
                  instances=[]
                  target=[]
              i=i+1
          else: 
             instances.append(np.array(kmer(lines[i].split(',')[0])))  # remove tabs  
             target.append(float(lines[i].split(',')[1].strip('\n'))) # save labels in list   \
             i=i+1
   # print i         
    return protein,labels    
    
if __name__ == '__main__':
    
    
    bags_ds1=create_bags('S1.txt')
    bags_ds2=create_bags('S2.txt')
    bags=bags_ds1+bags_ds2

    classifier=linclass(epochs=5000, Lambda=0.0001)

    Folds=create_folds(bags,5)
    prob=[]
    target=[]
    a=[]
    for i in range(5):
        train=[bags[Folds[i].train_bags[ii]] for ii in range(len(Folds[i].train_bags))]
        test=[bags[Folds[i].test_bags[ii]] for ii in range(len(Folds[i].test_bags))]
        classifier.train(train,merge=True)
        score=classifier.test(test)
        labels=[]
###        
        labels+=[b.label for b in test]     
        prob.append(score)
        target.append(labels)
        

    sco=[]
    lab=[]    
    for i in range(5):
      for k in range(len(prob[i])):
          sco.append(prob[i][k])
          lab.append(target[i][k])

    fpr, tpr, thresholds = roc_curve(lab,sco) # plotting ROC 
    a.append(auc(fpr,tpr)) 
        
    plt.plot(fpr,tpr,marker='.')
    plt.xlabel('fpr')
    plt.ylabel('tpr')
    plt.grid() 
    print np.mean(np.array(a))     
#    precision, recall, _ = precision_recall_curve(lab,sco)
#    av= average_precision_score(lab, sco)
#    print av
#    plt.plot(precision, recall, marker='o')
#plt.legend(['llclasss auc=0.886','linclass auc=0.865'],loc=4)
   
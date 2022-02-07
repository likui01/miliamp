# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 01:36:58 2016

@author: likui
"""

from CDhitMIL import *
from amyloids import *
from customized import *
from sklearn.externals import joblib
def linear_svm(): # testing on DS2, training on DS1
    amino, target, amino1, target1=read_data()
    Xtr=Compute_features(amino)
    scr=svm.SVC(C=0.1,kernel='linear', class_weight='auto')         
    scr.fit(Xtr, target)
    dd=[] # calculated targets 
    
    prob=[]
    for kk in range(len(amino1)): # use sliding window to calculate  the score 
        ww=[]#ds=np.array(3)# to store the score + labels 
        for yy in range(len(amino1[kk])-6):
            aa=kmer(amino1[kk][yy:yy+6])
            ww.append(np.hstack([scr.decision_function(aa)[0],scr.predict(aa)])) # using SVM to predict  the label 
        ww=np.array(ww)
        ind=np.argmax(ww[:,0])  
        #seq=kmer(amino1[kk][ind:ind+6])
        #prob.append(clf.predict_proba(seq)[0][1])
        prob.append(ww[ind,0])
        dd.append(ww[ind,1])
    return target1, prob   

def Linear_MIL(): # testing on DS2 and training on DS1+DS2
    bags_ds1=create_bags('S1.txt')
    bags_ds2=create_bags('S2.txt')
#    bag3=create_bags('S33win.txt')
    bags=bags_ds1+bags_ds2

    classifier=linclass(epochs=5000, Lambda=0.0001)

    Folds=create_folds(bags,5)
    prob=[]
    target=[]
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
    return lab,sco
def Customized_classifier(): #testing on DS2, training on DS1+DS2+Muatants-aggrescan 
       test1=create_bags('S1.txt')
       test2=create_bags('S2.txt')
#       train1,folds=create_mutants()
       bag_aggre, bag_other=sep_aggre()
       bags=test1+test2

       classifier=linclassrank(epochs=5000,Lambda=0.001)#0.0001
#   

       Folds=create_folds(bags,5)
       
       prob=[]
       target=[]
       
       for i in range(5):
                train=[bags[Folds[i].train_bags[ii]] for ii in range(len(Folds[i].train_bags))]
                test=[bags[Folds[i].test_bags[ii]] for ii in range(len(Folds[i].test_bags))]
#                classifier=linclassrank(epochs=5000,Lambda=1000)#0.0001
                classifier.train(train,bag_other)
                score=classifier.test(test)
                labels=[]
        ###        
                labels+=[b.label for b in test]     
                prob.append(score)
                target.append(labels)
       joblib.dump(classifier,'amyloid_pred.pkl')
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
    lab1,prob=linear_svm() # simple SVM 
    lab2,sco1=Linear_MIL() # MIL
    lab3,sco2=Customized_classifier() #linclassrank classifier
    lab4,sco3=readwaltz()
    sco4=readscore('mettamyscore.txt')
    sco5=readscore('agreescanmax.txt')
    sco6=readscore('appnnscore.txt')
    aa=[]
    for i in [(lab1,prob,'.','-'),(lab2,sco1,'o','-'),(lab3,sco2,"3",'-.'),(lab1,sco4,'>','-'),(lab1,sco5,'v',':'),(lab1,sco6,'2','--'),(lab4,sco3,'+',':')]:
        
        fpr, tpr, thresholds = roc_curve(i[0],i[1]) # plotting ROC 
        plt.plot(fpr,tpr, marker=i[2],linestyle=i[3])
        a=auc(fpr,tpr)
        aa.append(a)

    plt.xlabel('Fpr')
    plt.ylabel('Tpr')         
    plt.xlim([-0.1,0.5])

    plt.legend(['Linear SVM:'+str(round(aa[0],3)*100),'MIL:'+str(round(aa[1],3)*100),'MIL-Rank:'+str(round(aa[2],3)*100),'MetAmyl:'+str(round(aa[3],3)*100),'Aggrescan:'+str(round(aa[4],3)*100),'APPNN:'+str(round(aa[5],3)*100),'Waltz:'+str(round(aa[6],3)*100)],loc=4)
    plt.grid()







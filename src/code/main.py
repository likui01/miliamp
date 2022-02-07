# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 21:46:33 2018

@author: farzeen
"""


from CDhitMIL import *
from amyloids import *
from customized import *
import numpy as np
import pickle 
from sklearn.externals import joblib
from classifiers import *
import numpy as np
from scipy.io import *
from bag import *
from llc import *
from cv import *
from classifiers import *
from roc import*
#from roc import *
from temp2 import *
def test_bags_mutants(LIST):   
    mutants=[]
    original=[]
    target=[]
    idds=[]
    target=None
    for item in LIST:
         for i in item[1::]:
              
              if i[0]+3>len(item[0]):
                d=item[0][i[0]-6:i[0]]
                sc = copy.copy(d)
                original.append(sc)
                d=d.strip(i[1])
           # b=i[2]
                d+=i[2]
              elif i[0]-3<0:
                 d=item[0][i[0]:i[0]+6]
                 sc = copy.copy(d)
                 original.append(sc)
                 d=i[2]+d[1::]   
              else:
                
                d=item[0][i[0]-3:i[0]+3]
               
        
                sc = copy.copy(d)
                original.append(sc)
                d= d[:2]+i[2]+d[3:]
              mutants.append(d)
              #target.append(i[3])
              idds.append(str(i[1]+str(i[0])+str(i[2])))
   # target = 2 * np.array(target) - 1  
    #if len(seq[0])==6:
             
    feat1=[]
    for i in range(len(mutants)):
         aa=np.empty((1,20))
         ab=np.hstack((np.array(kmer(mutants[i]))-np.array(kmer(original[i]))))
#         print ab
         aa=np.matrix(ab)
         #print np.shape(aa.T)
         feat1.append(aa)
  
    bags=[Bag() for i in range(len(feat1))]
    for b in range(len(bags)):
        
        bags[b].instances=np.array(feat1[b])
        bags[b].label=target
        bags[b].peta=[1.0]
        bags[b].id=idds[b]
    return bags

def test_hotspot_bags(sequence):    
    target = None
    if len(sequence)==6:
         feat1=[]
         aa=np.empty((1,20))
         ab=np.hstack((kmer(sequence)))
         aa=np.matrix(ab)
         feat1.append(aa)
         bags=Bag()
         bags.instances=np.array(feat1[0])
         bags.label=target
         bags.peta=[1.0]
             
             
    else:
        feat1=[]
        for yy in range(len(sequence)-6):
                aa=np.empty((1,20))
                aa=np.array(np.hstack((kmer(sequence[yy:yy+6]))))
               # print aa
                ab=np.matrix(aa)
                feat1.append(np.array(ab))      
        bags=[Bag() for i in range(len(feat1))]
        for b in range(len(bags)):        
            bags[b].instances=np.array(feat1[b])
            bags[b].label=target
            bags[b].peta=[1.0]        
    return bags 
def create_bags(sequence):
 
    seq=[sequence]     
    target = None
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
            bags[b].label=target
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
            bags[b].label=target
            bags[b].peta=[1.0]*(np.shape(feat1[b])[0])        
    return bags    
if __name__ == '__main__':
    input_seq='STTAVASTAVV'
    hot_bag=test_hotspot_bags(input_seq)
    hot_clf=joblib.load('hot_spot_MIL_rank.pkl')
    results=hot_clf.test(hot_bag)
    amy_bag=create_bags(input_seq)
    amy_clf=joblib.load('amyloid_pred.pkl')
    res=amy_clf.test(amy_bag)
    amylinRAT='KCNTATCATQRLANFLVRSSNNLGPVLPPTNVGSNTY'
    LIST=[[(amylinRAT),(18,'R','H',1),(23,'L','F',1),(26,'V','I',1)]]
    mut=test_bags_mutants(LIST)
    mut_clf=joblib.load('mute_pred.pkl')
    re=mut_clf.test(mut)

























# -*- coding: utf-8 -*-
"""
Created on Thu Sep 08 10:57:59 2016

@author: likui
"""

from scipy.sparse import *
import random
import numpy as np
from cv import *
from bag import Bag
from classifiers import*

def hloss_fun(y, scores,**kwargs):
    
    return np.max(np.vstack((1-y*scores.flatten(),np.zeros(len(scores)))),axis = 0)

    
def hloss_dfun(x, y, z, **kwargs):
    if y*z<1:
        return -y*x
    else:
        return 0


                
class linclassrank(ClassifierBase):
    def train(self, bags1,bag2, **kwargs):
        ## bag1~MIL BAG
        ### bag2~ delta features
       if 'Beta' in kwargs:
                Beta=kwargs['Beta']
       else:          #print Beta
         Beta=0.0
       print Beta
       sizx=np.shape(bags1[0].instances)[1]
       w=np.zeros(sizx)
       
     #  w=w[np.newaxis,:]
      # print w
       epochs=self.epochs
       lamb=self.Lambda
       bias=self.bias
       pos,neg=separate_bags(bags1)
       pos_bags1=[]
       neg_bags1=[]
       for ii in pos:
           pos_bags1.append(bags1[ii])
       for ij in neg:
            neg_bags1.append(bags1[ij])
            
       pos1,neg1=separate_bags(bag2)
       pos_bag2=[]
       neg_bag2=[]
       for ki in pos1:
           pos_bag2.append(bag2[ki])
       for jk in neg1:
           neg_bag2.append(bag2[jk])
       T=(len(pos_bag2)+len(neg_bag2))*epochs 
       
#       import pdb; pdb.set_trace()
       for t in range(1,T):
           eta=1.0/(lamb*t)
#           if t%epochs==0:
#               random.shuffle(pos_bags1)
#               random.shuffle(pos_bag2)
#               random.shuffle(neg_bags1)
#               random.shuffle(neg_bag2)
           if(t%2)==0:
             bag1_sel=pos_bags1[t/2%len(pos_bags1)]
             bag2_sel=pos_bag2[t/2%len(pos_bag2)]
#             print t/2%len(pos_bag2)
           else:
               bag1_sel=neg_bags1[int((((t-1.0)/2.0)+1)%len(neg_bags1))]
               bag2_sel=neg_bag2[int((((t-1.0)/2.0)+1)%len(neg_bag2))]
#               print int((((t-1.0)/2.0)+1)%len(neg_bag2))
               
           X=bag1_sel.instances
           scores=X.dot(w.T)+bias
           indmax=np.argmax(scores)
           fxi=scores[indmax]
           xi=bag1_sel.instances[indmax]
           y_I=bag1_sel.label
           
           w=(1-eta*lamb)*w
#           bias=(1-eta*lamb)*bias
           #print y_I*fxi
           if(y_I*fxi<1):
#               print "hi"
               w=w+(1-Beta)*(eta*y_I*xi)
               bias=bias+(1-Beta)*(eta*y_I)
           ########################################################
           
           
           dX=bag2_sel.instances
           value=dX.dot(w.T)
           y_j=bag2_sel.label  
           
           if(y_j*value<1):
#               print "h"
               w=w+Beta*(eta*y_j*dX)
               
               
           
       self.w=w
       self.bias=bias
#       return self.w
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
       
           
           
           
           
        
            
            
       
       
       
       



















                
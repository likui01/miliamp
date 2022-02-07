# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 13:07:56 2016

@author: likui
"""
import numpy as np
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
from roc import*
def seg(filename):
    file = open(filename,"r")
    lines1 = file.readlines()
    file.close()
    amino1=[]
    target=[]
    for i in range(1,69,2):
        amino1.append(lines1[i].strip('\n'))
    for i in range(34):    
        target.append(np.zeros((1,len(amino1[i]))))
     
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
       
    labels=[]   
    proteins=[]
    for i in range(len(amino1)):
        
        qw=[]
        w1=[]
        for jj in range(len(amino1[i])-6):
        
            
            aa=amino1[i][jj:jj+6] #+kmer(seq[i][yy:yy+7])
            qw.append(aa)
            bb=sum(target[i][0,jj:jj+6])
            w1.append(bb)
        indexes = [ii for ii,x in enumerate(w1) if x == 5]   
    #    w1=np.array(w1)
        for ik in indexes:
          for ij in range(-8,0):  
              if (ik-ij)<0:
                  pass
              elif(ik--ij)<0:
                  pass
              else:
                  w1[ik--ij]=1.0
              if(ik+ij)>=len(w1):
                   pass
              elif(ik-ij)>=len(w1):
                  pass
              else:
                   w1[ik-ij]=1.0
                                
                    
        for iii in range(len(w1)):
            if w1[iii]>1.0:
                w1[iii]=1.0
        proteins.append(qw)   
        labels.append(w1)
#    import pdb;pdb.set_trace()    
    return labels, proteins
def read(filename):
    file=open(filename,'r')
    lines=file.readlines()
    file.close()
    s=[]
    
    for i in range(len(lines)):
        s.append(float(lines[i].strip('\n')))
     
    
    return s 
def compute(score,labels):
  sc=[]
  for i in range(33):
        aq=[]
        aw=[]
        ww=[]
        for jj in range(len(labels[i])):
            aa=score[i+jj]
            ww.append(aa)
            if labels[i][jj]==1:
                aq.append(jj)
                aw.append(aa)
            
            else:
             pass
        
        for ie in (range(len(aq))):
           ww[aq[ie]]=max(aw)   
        sc.append(ww)            
   
       
  sco=[]
  lab=[]    
  for i in range(33):
      for k in range(len(sc[i])):
          sco.append(sc[i][k])
          lab.append(labels[i][k])
    
  return sco,lab    
if __name__ == '__main__':
  labels,pro=seg('S333.txt')  
  aa=[]
  for i in [('linerasvm.txt',':','o'),('linearmil24.txt','--','2'),('ampmil.txt','-.','*'),('me.txt','-','>'),('aggre.txt',':','v'),('mett2.txt','--','4')]:
      s=read(i[0])    
      sco,lab=compute(s,labels)
              
      fpr, tpr, thresholds = roc_curve(lab,sco) # plotting ROC 
      a=auc(fpr,tpr)
      aa.append(a)
#      print a
      plt.plot(fpr,tpr, marker=i[2],linestyle=i[1])
#  plt.figure()      
  plt.xlabel('Fpr')  
  plt.ylabel('Tpr')
  plt.xlim([-0.05,0.4])   
  plt.grid() 
  plt.legend(['Linear SVM:'+str(round(aa[0],3)*100),'MIL:'+str(round(aa[1],3)*100),'MIL-Rank:'+str(round(aa[2],3)*100),'MetAmyl:'+str(round(aa[3],3)*100),'Aggrescan:'+str(round(aa[4],3)*100),'APPNN:'+str(round(aa[5],3)*100)],loc=4)
    
  #plt.legend([', auc='+str(round(a[0],3)),'APPNN, auc='+str(round(a[1],3)),'Aggrescan, auc='+str(round(a[2],3)),'MetAmyl, auc='+str(round(a[3],3)),'Linear MIL, train ds1+ds2, auc='+str(round(a[4],3)),'LLCMIL, train ds1+ds2, auc='+str(round(a[5],3)),'LLCMIL5cv, auc='+str(round(a[6],3))],loc=4)
           
       
       
       
       
   
#   
##   
#with open('S33w.fasta','w')as f:
#    for i in range(len(pro)-1):
#        for j in range(len(pro[i])): 
#             f.write('>tro|'+str(i)+"|"+str(labels[i][j]))
#             f.write('\n')
#             f.write(str(pro[i][j]))
#            
#             f.write("\n")
#        
#f.close()        
##   
##   
#   
#   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
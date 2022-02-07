# -*- coding: utf-8 -*-
"""
Created on Fri Sep 02 10:24:11 2016

@author: likui
"""
from classifiers import *
import numpy as np
from scipy.io import *
from bag import *
from llc import *
from cv import *
#from classifiersrc import *
#from classifiers import *
from mute import *
from roc import*
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
import matplotlib.pyplot as plt #importing plotting module
#from roc import *
import copy
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
             
def sep_aggre():
    
    amyloidab='DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA'
    amyseg1='VHHQKLVFFAEDVGS'
    amyseg2='KKLVFFAED'
    alpha='MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA'
    Acp='STAQSLKSVDYEVFGRVQGVSFRMYTEDEARKIGVVGWVKNTSKGTVTGQVQGPEDKVNSMKSWLSKVGSPSSRIDRTNFSNEKTISKLEYSNFSVRY'
    TTR='GPTGTGESKCPLMVKVLDAVRGSPAINVAVHVFRKAADDTWEPFASGKTSESGELHGLTTEEEFVEGIYKVEIDTKSYWKALGISPFHEHAEVVFTANDSGPRRYTIAALLSPYSYSTTAVVTNPKE'
    tau='MAEPRQEFEVMEDHAGTYGLGDRKDQGGYTMHQDQEGDTDAGLKESPLQTPTEDGSEEPGSETSDAKSTPTAEDVTAPLVDEGAPGKQAAAQPHTEIPEGTTAEEAGIGDTPSLEDEAAGHVTQARMVSKSKDGTGSDDKKAKGADGKTKIATPRGAAPPGQKGQANATRIPAKTPPAPKTPPSSGEPPKSGDRSGYSSPGSPGTPGSRSRTPSLPTPPTREPKKVAVVRTPPKSPSSAKSRLQTAPVPMPDLKNVKSKIGSTENLKHQPGGGKVQIINKKLDLSNVQSKCGSKDNIKHVPGGGSVQIVYKPVDLSKVTSKCGSLGNIHHKPGGGQVEVKSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETHKLTFRENAKAKTDHGAEIVYKSPVVSGDTSPRHLSNVSSTGSIDMVDSPQLATLADEVSASLAKQGL'    
    amylinhuman='KCNTATCATQRLANFLVHSSNNFGAILSSTNVGSNTY' 
    SOD1='ATKAVCVLKGDGPVQGIINFEQKESNGPVKVWGSIKGLTEGLHGFHVHEFGDNTAGCTSAGPHFNPLSRKHGGPKDEERHVGDLGNVTADKDGVADVSIEDSVISLSGDHCIIGRTLVVHEKADDLGKGGNEESTKTGNAGSRLACGVIGIAQ'
    PrP='MANLGCWMLVLFVATWSDLGLCKKRPKPGGWNTGGSRYPGQGSPGGNRYPPQGGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQGGGTHSQWNKPSKPKTNMKHMAGAAAAGAVVGGLGGYMLGSAMSRPIIHFGSDYEDRYYRENMHRYPNQVYYRPMDEYSNQNNFVHDCVNITIKQHTVTTTTKGENFTETDVKMMERVVEQMCITQYERESQAYYQRGSSMVLFSSPPVILLISFLIFLIVG'
    amylinRAT='KCNTATCATQRLANFLVRSSNNLGPVLPPTNVGSNTY' 
    LIST=[]#aggrescan mutants 
    LIST.append([(amyloidab),(21,'A','G',-1),(22,'E','K',1),(22,'E','G',1),(22,'E','Q',1),(19,'F','P',-1),(19,'F','T',-1),(23,'D','N',1),(19,"F",'D',-1),(31,'I','L',-1),(32,'I','L',-1),(41,'I','G',-1),(41,'I','A',-1),(41,'I','L',-1),(42,'A','G',-1),(42,'A','V',1)])
    LIST.append([(alpha),(46,'E','K',1),(53,'A','T',1),(76,'A','E',-1),(76,'A','R',-1)])
    LIST.append([(tau),(5,'R','L',1),(406,'R','W',1),(272,'G','V',1),(320,'S','F',1),(301,'P','L',1)])
    LIST.append([(amylinhuman),(22,'N','A',1),(23,'F','A',-1),(24,'G','A',1),(26,'I','A',-1),(27,'L','A',-1),(20,'S','G',1)])
    #LIST.append([(PrP),(111,'H','A',1),(111,'H','K',-1),(117,'A','V',1),(210,'V','I',1)])
    LIST.append([(amylinRAT),(18,'R','H',1),(23,'L','F',1),(26,'V','I',1)])
    LIST2=[]# other mutants 
    LIST2.append([(amyseg1),(2,'H','P',-1),(3,'H','P',-1),(4,'Q','P',-1),(5,'K','P',-1),(8,'F','P',-1),(11,'E','P',-1),(12,'D','P',-1),(13,'V','P',-1),(14,'G','P',-1)]) 
    LIST2.append([(amyseg2),(3,'L','P',-1),(4,'V','P',-1),(5,'F','P',-1),(6,'F','P',-1),(7,'A','P',-1)])
    LIST2.append([(Acp),(5,'S','T',-1),(9,'V','A',1),(11,'Y','F',1),(13,'V','A',-1),(17,'V','A',-1),(20,'V','A',-1),(22,'F','L',-1),(25,'Y','A',-1),(29,'E','D',1),(30,'A','G',1),(33,'I','V',1),(34,'G','A',-1),(36,'V','A',-1),(39,'V','A',-1),(47,'V','A',-1),(51,'V','A',-1),(54,'P','A',1),(61,'M','A',1),(64,'W','F',-1),(65,'L','V',-1),(76,'K','A',-1),(71,'P','A',1),(75,'I','V',-1),(78,'S','T',-1),(83,'E','D',-1),(86,'I','V',1),(87,'S','T',1),(89,'L','A',-1),(91,'Y','Q',-1),(92,'S','T',1),(94,'F','L',-1),(98,'Y','Q',-1),(8,'S','H',1),(21,'S','R',1),(23,'R','Q',1),(29,'E','K',1),(29,'E','R',-1),(43,'S','E',1),(55,'E','Q',-1),(77,'R','E',1),(88,'K','N',1),(88,'K','Q',1),(90,'E','H',1),(92,'S','H',1),(97,'R','E',1),(97,'R','Q',1)]) 
    LIST2.append([(tau),(5,'R','H',-1),(257,'K','T',-1),(266,'L','V',1),(279,'N','K',1),(296,'N','H',-1),(301,'P','S',1),(305,'S','N',-1),(389,'G','R',-1),(337,'V','M',1),(342,'E','V',-1),(369,'K','I',1)])
    #LIST2.append([(TTR),(30,'V','M',1),(119,'T','M',-1),(55,'L','P',1)])
    LIST2.append([(SOD1),(4,'A','V',1),(4,'A','T',1),(93,'G','A',1),(93,'G','D',1),(93,'G','V',1),(84,'L','V',1),(90,'A','D',1),(124,'D','A',1),(14,'V','A',1),(14,'V','M',1),(21,'E','K',1),(41,'G','S',1),(41,'G','D',1),(100,'E','G',1),(100,'E','K',1),(139,'N','K',1),(93,'G','C',1),(43,'H','R',1),(101,'D','N',1),(101,'D','G',1),(144,'L','F',1),(144,'L','S',1),(148,'V','G',1),(148,'V','I',1)])
    bag_aggre,f=create_bags_mutants(LIST,make_fold=False)    
    bag_other,f=create_bags_mutants(LIST2,make_fold=False)
    return bag_aggre,bag_other
    
    
        
      
      
def create_mutants():    
    amyloidab='DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA'
    amyseg1='VHHQKLVFFAEDVGS'
    amyseg2='KKLVFFAED'
    alpha='MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA'
    Acp='STAQSLKSVDYEVFGRVQGVSFRMYTEDEARKIGVVGWVKNTSKGTVTGQVQGPEDKVNSMKSWLSKVGSPSSRIDRTNFSNEKTISKLEYSNFSVRY'
    TTR='GPTGTGESKCPLMVKVLDAVRGSPAINVAVHVFRKAADDTWEPFASGKTSESGELHGLTTEEEFVEGIYKVEIDTKSYWKALGISPFHEHAEVVFTANDSGPRRYTIAALLSPYSYSTTAVVTNPKE'
    tau='MAEPRQEFEVMEDHAGTYGLGDRKDQGGYTMHQDQEGDTDAGLKESPLQTPTEDGSEEPGSETSDAKSTPTAEDVTAPLVDEGAPGKQAAAQPHTEIPEGTTAEEAGIGDTPSLEDEAAGHVTQARMVSKSKDGTGSDDKKAKGADGKTKIATPRGAAPPGQKGQANATRIPAKTPPAPKTPPSSGEPPKSGDRSGYSSPGSPGTPGSRSRTPSLPTPPTREPKKVAVVRTPPKSPSSAKSRLQTAPVPMPDLKNVKSKIGSTENLKHQPGGGKVQIINKKLDLSNVQSKCGSKDNIKHVPGGGSVQIVYKPVDLSKVTSKCGSLGNIHHKPGGGQVEVKSEKLDFKDRVQSKIGSLDNITHVPGGGNKKIETHKLTFRENAKAKTDHGAEIVYKSPVVSGDTSPRHLSNVSSTGSIDMVDSPQLATLADEVSASLAKQGL'    
    amylinhuman='KCNTATCATQRLANFLVHSSNNFGAILSSTNVGSNTY' 
    SOD1='ATKAVCVLKGDGPVQGIINFEQKESNGPVKVWGSIKGLTEGLHGFHVHEFGDNTAGCTSAGPHFNPLSRKHGGPKDEERHVGDLGNVTADKDGVADVSIEDSVISLSGDHCIIGRTLVVHEKADDLGKGGNEESTKTGNAGSRLACGVIGIAQ'
    PrP='MANLGCWMLVLFVATWSDLGLCKKRPKPGGWNTGGSRYPGQGSPGGNRYPPQGGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQPHGGGWGQGGGTHSQWNKPSKPKTNMKHMAGAAAAGAVVGGLGGYMLGSAMSRPIIHFGSDYEDRYYRENMHRYPNQVYYRPMDEYSNQNNFVHDCVNITIKQHTVTTTTKGENFTETDVKMMERVVEQMCITQYERESQAYYQRGSSMVLFSSPPVILLISFLIFLIVG'
    amylinRAT='KCNTATCATQRLANFLVRSSNNLGPVLPPTNVGSNTY' 
    LIST=[]
    LIST.append([(amyloidab),(21,'A','G',-1),(22,'E','K',1),(22,'E','G',1),(22,'E','Q',1),(19,'F','P',-1),(19,'F','T',-1),(23,'D','N',1),(19,"F",'D',-1),(31,'I','L',-1),(32,'I','L',-1),(41,'I','G',-1),(41,'I','A',-1),(41,'I','L',-1),(42,'A','G',-1),(42,'A','V',1)])
    LIST.append([(amyseg1),(2,'H','P',-1),(3,'H','P',-1),(4,'Q','P',-1),(5,'K','P',-1),(8,'F','P',-1),(11,'E','P',-1),(12,'D','P',-1),(13,'V','P',-1),(14,'G','P',-1)]) 
    LIST.append([(amyseg2),(3,'L','P',-1),(4,'V','P',-1),(5,'F','P',-1),(6,'F','P',-1),(7,'A','P',-1)])
    LIST.append([(Acp),(5,'S','T',-1),(9,'V','A',1),(11,'Y','F',1),(13,'V','A',-1),(17,'V','A',-1),(20,'V','A',-1),(22,'F','L',-1),(25,'Y','A',-1),(29,'E','D',1),(30,'A','G',1),(33,'I','V',1),(34,'G','A',-1),(36,'V','A',-1),(39,'V','A',-1),(47,'V','A',-1),(51,'V','A',-1),(54,'P','A',1),(61,'M','A',1),(64,'W','F',-1),(65,'L','V',-1),(76,'K','A',-1),(71,'P','A',1),(75,'I','V',-1),(78,'S','T',-1),(83,'E','D',-1),(86,'I','V',1),(87,'S','T',1),(89,'L','A',-1),(91,'Y','Q',-1),(92,'S','T',1),(94,'F','L',-1),(98,'Y','Q',-1),(8,'S','H',1),(21,'S','R',1),(23,'R','Q',1),(29,'E','K',1),(29,'E','R',-1),(43,'S','E',1),(55,'E','Q',-1),(77,'R','E',1),(88,'K','N',1),(88,'K','Q',1),(90,'E','H',1),(92,'S','H',1),(97,'R','E',1),(97,'R','Q',1)]) 
    LIST.append([(alpha),(30,"A",'P',-1),(46,'E','K',1),(53,'A','T',1),(76,'A','E',-1),(76,'A','R',-1)])
    LIST.append([(tau),(5,'R','L',1),(406,'R','W',1),(272,'G','V',1),(310,'Y','W',1),(320,'S','F',1),(301,'P','L',1),(5,'R','H',-1),(257,'K','T',-1),(266,'L','V',1),(279,'N','K',1),(296,'N','H',-1),(301,'P','S',1),(305,'S','N',-1),(389,'G','R',-1),(337,'V','M',1),(342,'E','V',-1),(369,'K','I',1)])
    LIST.append([(TTR),(30,'V','M',1),(119,'T','M',-1),(55,'L','P',1)])
    LIST.append([(amylinhuman),(22,'N','A',1),(23,'F','A',-1),(24,'G','A',1),(26,'I','A',-1),(27,'L','A',-1),(20,'S','G',1)])
    LIST.append([(SOD1),(4,'A','V',1),(4,'A','T',1),(93,'G','A',1),(93,'G','D',1),(93,'G','V',1),(84,'L','V',1),(90,'A','D',1),(124,'D','A',1),(14,'V','A',1),(14,'V','M',1),(21,'E','K',1),(41,'G','S',1),(41,'G','D',1),(100,'E','G',1),(100,'E','K',1),(139,'N','K',1),(93,'G','C',1),(43,'H','R',1),(101,'D','N',1),(101,'D','G',1),(144,'L','F',1),(144,'L','S',1),(148,'V','G',1),(148,'V','I',1)])
    LIST.append([(PrP),(111,'H','A',1),(111,'H','K',-1),(117,'A','V',1),(210,'V','I',1)])
    LIST.append([(amylinRAT),(18,'R','H',1),(23,'L','F',1),(26,'V','I',1)])
    bags,folds=create_bags_mutants(LIST,make_fold=True)
    return bags,folds
def create_bags_mutants(LIST,make_fold=True):   
    mutants=[]
    original=[]
    target=[]
    idds=[]
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
              target.append(i[3])
              idds.append(str(i[1]+str(i[0])+str(i[2])))
   # target = 2 * np.array(target) - 1  
    #if len(seq[0])==6:
    folds=[]          
    if make_fold==True:
         folds=[fold() for _ in range(9)]   
         aaa=len(LIST[0])+len(LIST[1])+len(LIST[2])-3
         index=aaa
         ind=[]
         ind.append(aaa)
         for iii in range(3,11):
              
              index=index+len(LIST[iii])-1
              ind.append(index)
         folds[0].test_bags=range(aaa)
         folds[0].train_bags=range(aaa,len(mutants))
         for i in range(1,9):
             folds[i].test_bags=range(ind[i-1],ind[i])
             folds[i].train_bags=range(ind[i-1])+range(ind[i],len(mutants))
            
             
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
        bags[b].label=target[b]
        bags[b].peta=[1.0]
        bags[b].id=idds[b]
    return bags, folds
if __name__ == '__main__':  
       test1=create_bags('S1.txt')
       test2=create_bags('S2.txt')
       train1,folds=create_mutants()
       bag_aggre, bag_other=sep_aggre()
       bag3=create_bags('S33win.txt')
       bags=test1+test2
#       for i in range(10):
#       compute_gammas(BAGS, C=40, k=40, beta=1)
    #    compute_gammas(bags_ds3,C=50,k=50,beta=10.0)
#       classifier=llclass(epochs=5000, Lambda=0.001)
#       classifier=linclass(epochs=5000, Lambda=0.001)
    #    tes,target=readtestinstance('33window.txt')
#       classifier=linclass_rank(epochs=5000,Lambda=0.1)#,merge=True)
#       classifier=llclass_rank(epochs=5000,Lambda=0.0001)#merge=True)
       classifier=linclassrank(epochs=5000,Lambda=0.00001)#0.0001
#       test=BAGS[0:len(test1)]
#       train=BAGS[len(test1)::]
########################################################################################################3
       prob=[]
       target=[]
##       test=[]
##       train=[]
###       for i in [0,2,3,5]:# range(9):
###           
###            test+=[train1[folds[i].test_bags[ii]] for ii in range(len(folds[i].test_bags))]
###       for jj in [1,4,6,7,8]:    
###            train+=[train1[folds[jj].test_bags[ii]] for ii in range(len(folds[jj].test_bags))]
###                classifier=linclassrank(epochs=5000,Lambda=1000)#0.0001
##       for i in range(9):
##           test=[train1[folds[i].test_bags[ii]] for ii in range(len(folds[i].test_bags))]
##           
##           train=[train1[folds[i].train_bags[ii]] for ii in range(len(folds[i].train_bags))]
#        
       classifier.train(bags,bag_other)
       score=classifier.test(bag_aggre)
       labels=[]
       ids=[]
#    ###        
       labels+=[b.label for b in bag_aggre]    
       ids+=[b.id for b in bag_aggre]
       prob.append(score)
       target.append(labels)
#       sco=[]
#       lab=[]    
#       for i in range(9):
#          for k in range(len(prob[i])):
#              sco.append(prob[i][k])
#              lab.append(target[i][k])    
  ################################################################################################## 
#       Folds=create_folds(bags,5)
#       
#       prob=[]
#       target=[]
#       a=[]
#       weig=[]
#       for i in range(5):
#                train=[bags[Folds[i].train_bags[ii]] for ii in range(len(Folds[i].train_bags))]
#                test=[bags[Folds[i].test_bags[ii]] for ii in range(len(Folds[i].test_bags))]
##                classifier=linclassrank(epochs=5000,Lambda=1000)#0.0001
#                classifier.train(train,bag_other)
#                score=classifier.test(test)
##                prob.append(score)
#                labels=[]
#        ###        
#                labels+=[b.label for b in test]     
#                prob.append(score)
#                target.append(labels)
#                weig.append(classifier.w)
##       classifier.train(bags,train1)
##       score=classifier.test(train1)
##       labels=[]
##    ###        
##       labels+=[b.label for b in train]     
#       sco=[]
#       lab=[]    
#       for i in range(5):
#          for k in range(len(prob[i])):
#              sco.append(prob[i][k])
#              lab.append(target[i][k])    
##      
##    
       fpr, tpr, thresholds = roc_curve(labels,score) # plotting ROC 
       a=auc(fpr,tpr)
#        
#       plt.plot(fpr,tpr,marker='.')
#       plt.xlabel('fpr')
#       plt.ylabel('tpr')
#       plt.grid() 
#       print np.mean(np.array(a))
       
       
#
#er=np.mean(np.array(weig),axis=0).T
er=(classifier.w.T)-np.mean(classifier.w)     
sd=np.array(['A', 'C', 'E', 'D', 'G', 'F', 'I', 'H', 'K', 'M', 'L', 'N', 'Q', 'P', 'S', 'R', 'T', 'W', 'V', 'Y'])
dr=np.array([-0.036,0.604,-1.412,-1.832,-0.535,1.754,1.822,-1.033,-0.931,0.910,1.380,-1.302,-1.23,-0.33,-1.240,-0.159,1.037, 1.594,1.159]).T             
aa=zip(dr,sd)
#ssd=np.argsort(sd)
#er=classifier.w.T
we=[]
for i in range(20):
    we.append(er[i])
    
#aw=np.sort(sd)
nd=range(20)
plt.bar(nd,er)
#plt.grid()
bar_width=0.35
nd=np.array(nd)+0.5
plt.xticks(nd,sd)

plt.xlabel('Amnio Acid')
plt.ylabel('Weights')  
plt.grid()  
plt.show()   
###       
#       with open('S33scor.txt','w')as f:
#            for i in range(len(score)):
#                f.write(str(score[i]))
#                f.write("\n")
#                
#       f.close()

##    
#all2 = zip(score,labels,ids)
#all2.sort(key=lambda score: score[0])        
#           
#aggre=[-16,15,29,5,-68,-63,16,-118,-15,-15,-62,-49,-12,-10,32,2,2,2,1,2,2,-1,-5,-3,9,17,11,21,-59,16,-61,-23,-106]
#tar=[-1,1,1,1,-1,-1,1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1,-1,-1,1,1,1,1,-1,1,-1,-1,1]      
##aa3=zip(aggre,tar)
##aa3.sort(key=lambda aggre: aggre[0])   
#fpr, tpr, thresholds = roc_curve(tar,aggre) # plotting ROC 
#a=auc(fpr,tpr)
#
#plt.plot(fpr,tpr,marker='o')
#   plt.xlabel('Fpr')
#   plt.ylabel('Tpr')
#   plt.grid() 
#   print np.mean(np.array(a))  
#
#f = open('linearmil.txt', 'w')
#for ii in range(len(tt)):
#    f.write(str(tt[ii]))
#    f.write('\n')
#f.close()  
#####  
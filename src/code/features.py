# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 03:43:03 2016

@author: likui
"""



import numpy as np
import re
import ast
from Bio.Data import IUPACData as ip
def AAindex1():
    file = open('aaindex1.txt',"r") #open data file
    lines = file.readlines() # readlines
    file.close()
    n=[]
    for i in range(len(lines)):
        if (lines[i][0]=='I'):
            n.append(lines[i+1])
            n.append(lines[i+2])
        else:
            pass#  
    dictt= {'A':[],
            'R':[],
            'N':[],
            'D':[],
            'C':[],
            'Q':[],
            'E':[],
            'G':[],
            'H':[],
            'I':[],
            'L':[],
            'K':[],
            'M':[],
            'F':[],
            'P':[],
            'S':[],
            'T':[],
            'W':[],
            'Y':[],
            'V':[],}
    # remove NA values        
    n1=[]               
    for i in range(len(n)):
        if n[i].find('NA')>=0:
            if (i % 2 == 0): #even 
                 pass
            else: #odd
              n1.pop((len(n1)-1))
        else:
             n1.append(n[i])
             
    a=[]
      
    for i in range(len(n1)):
            n1[i]=(n1[i].replace('\t ',',').replace('\n','').split(','))  
            d=re.sub("\s\s+" , " ", n1[i][0])      
            d=d.replace('','').split(',') 
            f=[ ", ".join(item.split(" ")) for item in d ]
            a.append([float(k) for k in f[0][1:].split(',')])
    ## normalize         
    dd=[]
    for i in range(0,len(a)/2,1):
        dd.append(a[i]+a[i+1])
        mim=min(dd[i])
        maxx=max(dd[i])
        d=maxx-mim
        dd[i]=[ (dd[i][ii]-mim)/d for ii in range(len(dd[i])) ]
        
        
        
    for i in range(len(dd)):     
            dictt['A'].append(dd[i][0])
            dictt['R'].append(dd[i][1])
            dictt['N'].append(dd[i][2])
            dictt['D'].append(dd[i][3])
            dictt['C'].append(dd[i][4])
            dictt['Q'].append(dd[i][5])
            dictt['E'].append(dd[i][6])
            dictt['G'].append(dd[i][7])
            dictt['H'].append(dd[i][8])
            dictt['I'].append(dd[i][9])
            dictt['L'].append(dd[i][10])
            dictt['K'].append(dd[i][11])
            dictt['M'].append(dd[i][12])
            dictt['F'].append(dd[i][13])
            dictt['P'].append(dd[i][14])
            dictt['S'].append(dd[i][15])
            dictt['T'].append(dd[i][16])
            dictt['W'].append(dd[i][17])
            dictt['Y'].append(dd[i][18])
            dictt['V'].append(dd[i][19])    
    return dictt


def APDbase():
    file = open('APDbase.txt',"r") #open data file
    lines = file.readlines() # readlines
    file.close()
    n=[]
    dictt= {'A':[],
            'R':[],
            'N':[],
            'D':[],
            'C':[],
            'Q':[],
            'E':[],
            'G':[],
            'H':[],
            'I':[],
            'L':[],
            'K':[],
            'M':[],
            'F':[],
            'P':[],
            'S':[],
            'T':[],
            'W':[],
            'Y':[],
            'V':[],}
    for i in range(len(lines)):
                        n.append(lines[i].replace('\t',',').replace('\n',',').split(','))  # remove tabs 
             
                
                        dictt['A'].append(ast.literal_eval(n[i][0][1:]))
                        dictt['R'].append(float(n[i][1]))
                        dictt['N'].append(float(n[i][2]))
                        dictt['D'].append(float(n[i][3]))
                        dictt['C'].append(float(n[i][4]))
                        dictt['Q'].append(float(n[i][5]))
                        dictt['E'].append(float(n[i][6]))
                        dictt['G'].append(float(n[i][7]))
                        dictt['H'].append(float(n[i][8]))
                        dictt['I'].append(float(n[i][9]))
                        dictt['L'].append(float(n[i][10]))
                        dictt['K'].append(float(n[i][11]))
                        dictt['M'].append(float(n[i][12]))
                        dictt['F'].append(float(n[i][13]))
                        dictt['P'].append(float(n[i][14]))
                        dictt['S'].append(float(n[i][15]))
                        dictt['T'].append(float(n[i][16]))
                        dictt['W'].append(float(n[i][17]))
                        dictt['Y'].append(float(n[i][18]))
                        dictt['V'].append(float(n[i][19][:-1]))   
    return dictt
                  
    
    
def encode(): 
    
    dictt= {'A':[ 1 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ],
            'R':[ 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ],
            'N':[ 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ],
            'D':[ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ],
            'C':[ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ],
            'Q':[ 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ],
            'E':[ 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ],
            'G':[ 0, 0, 0, 0, 0 ,0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ],
            'H':[ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 ,0 ,0, 0, 0, 0, 0, 0, 0, 0, 0 ],
            'I':[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ],
            'L':[ 0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ],
            'K':[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 ],
            'M':[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 ],
            'F':[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 ],
            'P':[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 ],
            'S':[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0, 0, 0, 0, 1, 0, 0 ,0, 0 ],
            'T':[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 ,0, 0, 0 ],
            'W':[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 ],
            'Y':[ 0 ,0 ,0 ,0 ,0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 ],
            'V':[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0, 0, 0, 0, 0, 0, 0, 1 ],}
            
    
   
    #print 'Amino Acid Counts: ', counts    
    return dictt   



def kmer(seq):         # function to count 1 mer 
    counts = {}
    for k in ip.protein_letters: counts[k] = 0
    for a in ip.protein_letters:
        if a in seq:
            counts[a] = np.float(seq.count(a))

    #print 'Amino Acid Counts: ', counts    
    return counts.values()   # return the number of counts  


def encode1(seq): 
    
    dictt= {'A':[ 1 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ],
            'R':[ 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ],
            'N':[ 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ],
            'D':[ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ],
            'C':[ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ],
            'Q':[ 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ],
            'E':[ 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ],
            'G':[ 0, 0, 0, 0, 0 ,0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ],
            'H':[ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 ,0 ,0, 0, 0, 0, 0, 0, 0, 0, 0 ],
            'I':[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ],
            'L':[ 0, 0, 0, 0, 0, 0, 0 ,0 ,0 ,0 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ],
            'K':[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 ],
            'M':[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 ],
            'F':[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0 ],
            'P':[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0 ],
            'S':[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0, 0, 0, 0, 1, 0, 0 ,0, 0 ],
            'T':[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 ,0, 0, 0 ],
            'W':[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0 ],
            'Y':[ 0 ,0 ,0 ,0 ,0 ,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0 ],
            'V':[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ,0, 0, 0, 0, 0, 0, 0, 1 ],}
            
    AA=[]
    for k in range(len(seq)):
        AA=AA+dictt[seq[k]]
   
    #print 'Amino Acid Counts: ', counts    
    return AA    




def encode2(seq): 
    
    dictt1= {'A':[ -1.56,-1.67,-0.97,-.27,-.93,-0.78,-.2,-0.08,0.21,-0.48 ],
            'R':[0.22,1.27,1.37,1.87,-1.70,0.46,0.92,-.39,.23,0.93  ],
            'N':[ 1.14,-0.07,-.12,0.81,0.18,0.37,-0.09,1.23,1.10,-1.73],
            'D':[ 0.58,-.22,-1.58,0.81,-0.92,0.15,-1.52,0.47,0.76,0.70],
            'C':[ 0.12,-.89,0.45,-1.05,-0.71,2.41,1.52,-.69,1.13,1.1 ],
            'Q':[ -.47,0.24,0.07,1.1,1.1,0.59,0.84,-0.71,-0.03,-2.33],
            'E':[-1.45,0.19,-1.61,1.17,-1.31,0.4,0.04,0.38,-.35,-0.12 ],
            'G':[ 1.46,-1.96,-0.23,-0.16,0.1,-.11,1.32,2.36,-1.66,0.46 ],
            'H':[-0.41,0.52,-0.28,0.28,1.61,1.01,-1.85,0.47,1.13,1.63 ],
            'I':[ -.73,-.16,1.79,-0.77,-0.54,0.03,-0.83,0.51,0.66,-1.78,],
            'L':[ -1.04,0,-0.24,-1.1,-.55,-2.05,0.96,-0.76,0.45,0.93],
            'K':[ -.34,0.82,-.23,1.7,1.54,-1.62,1.15,-.08,-0.48,0.6 ],
            'M':[ -1.4,0.18,-.42,-0.73,2,1.52,.26,0.11,-1.27,0.27],
            'F':[ -.21,0.98,-.36,-1.43,0.22,-0.81,0.67,1.1,1.71,-.44 ],
            'P':[ 2.06,-.33,-1.15,-0.75,0.88,-0.45,0.3,-2.3,0.74,-0.28],
            'S':[ 0.81,-1.08,0.16,0.42,-0.21,-0.43,-1.89,-1.15,-0.97,-.23 ],
            'T':[ .26,-0.7,1.21,0.63,-0.1,0.21,0.24,-1.15,-0.56,0.19 ],
            'W':[ 0.3,2.1,-0.72,-1.57,-1.16,0.57,-0.48,-0.40,-2.30,-.6],
            'Y':[ 1.38,1.48,0.8,-0.56,-0,-0.68,-0.31,1.03,-0.05,0.53 ],
            'V':[ -0.74,-.71,2.04,-0.4,0.5,-.81,-1.07,0.06,-.46,0.65 ],}
            
    AA1=[]
    for k in range(len(seq)):
        AA1=AA1+dictt1[seq[k]]
   
    #print 'Amino Acid Counts: ', counts    
    return AA1   


def weightmer(seq):
    dictt= {'A':-0.036,
            'R':-1.240,
            'N':-1.302,
            'D':-1.836,
            'C':0.604,
            'Q':-1.231,
            'E':-1.412,
            'G':-0.535,
            'H':-1.033,
            'I':1.822,
            'L':1.380,
            'K':-0.931,
            'M':0.910,
            'F':1.754,
            'P':-0.334,
            'S':-0.294,
            'T':-0.159,
            'W':1.037,
            'Y':1.159,
            'V':1.594,}
    dic= {'A':8.25,
            'R':5.53,
            'N':4.06,
            'D':5.45,
            'C':1.37,
            'Q':3.39,
            'E':6.75,
            'G':7.07,
            'H':2.27,
            'I':5.96,
            'L':9.66,
            'K':5.84,
            'M':2.42,
            'F':3.86,
            'P':4.70,
            'S':6.56,
            'T':5.34,
            'W':1.08,
            'Y':2.92,
            'V':6.87,}      
    counts = {}
    for k in ip.protein_letters: counts[k] = 0
    for a in ip.protein_letters:
        if a in seq:
            counts[a] = seq.count(a)*float(dictt[a])#/dic[a]

    #print 'Amino Acid Counts: ', counts    
    return sum(counts.values())   # return the number of counts          

#encode2('A')

      
#for i in range(544):
#  k=0  
#  for j in range(20):
#      if j<10:
#        ss[j][i]=a[i][j]
#      else:
#        ss[j][i]=a[i+1][k]
#        k=k+1
#    #return ss        
#    
#    
#    
#sd=AAindex1()    
# 
#    
#    
#    
#    
#    
##a[].index(NA)    
#x = np.array(range(20))
##y = np.array([20,21,22,23])
#plt.xticks(x, n)
#plt.plot(x, clf.coef_.T, marker='o')
#plt.show()
#    
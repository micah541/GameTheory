import numpy as np
import scipy.special as sp
from scipy.stats import nbinom


### doublespend


p = 0.3
eps = 0.2
q = 1- p 
k = 6
### charlie reorg

def CharlieRound(l,k, p): 
    state_vector1 = 1
    state_vector2 = 0
    while(1):
        if state_vector1 == -1:
            return(l)
        if state_vector1 == k:
            return(-state_vector2)
        if np.random.uniform()<p:
            state_vector2 += 1
            state_vector1 += -1 
        else : 
            state_vector1 += 1
    

def CharlieSim(l,k,p,N):
    v = 0 
    for i in range(N): 
        v = v+CharlieRound(l,k,p)
    return(np.float(v)/N)

### the break even table checks out!! 


### Now check the e function for other values

def e_round(l,k,p,x):
    state_vector1 = x
    state_vector2 = 0
    while(1):
        if state_vector1 == -1:
            return(l)
        if state_vector1 == k:
            return(-state_vector2)
        if np.random.uniform()<p:
            state_vector2 += 1
            state_vector1 += -1 
        else : 
            state_vector1 += 1
    

def eCharlieSim(l,k,p,x,N):
    v = 0 
    for i in range(N): 
        v = v+e_round(l,k,p,x)
    return(np.float(v)/N)
    
## what my constants say it should be 


def e_compute(l,k,p,x):
    m = p/(1-p)
    r = (1-p)/(1-2*p)
    M_0 =np.asarray( [[np.power(m,-1), 1], [np.power(m,k), 1]])
    v_1 = np.asarray([1,0])
    M_0i = np.linalg.inv(M_0)
    [c_1,c_2] = np.matmul(M_0i,v_1)
    M_11 = np.asarray( [[np.power(m,-1), m], [-k*np.power(m,k), -k*m]])
    [c_1tlambda, c_2tlambda]=np.matmul(M_0i,v_1)
    [c_1t, c_2t]=r*np.matmul(M_11,[c_1,1-c_2])
    







def CharlieRound(l,k, p): 
    state_vector1 = 1
    state_vector2 = 0
    while(1):
        if state_vector1 == -1:
            return(l)
        if state_vector1 == k:
            return(-state_vector2)
        if np.random.uniform()<p:
            state_vector2 += 1
            state_vector1 += -1 
        else : 
            state_vector1 += 1
    







def ratio_six_block_ds(p):
    m = p/(1-p)
    r = (1-p)/(1-2*p)
    k = 6
    M_0 =np.asarray( [[np.power(m,-1), 1], [np.power(m,k), 1]])
    v_1 = np.asarray([1,0])
    M_0i = np.linalg.inv(M_0)
    [c_1,c_2] = np.matmul(M_0i,v_1)
    M_11 = np.asarray( [[np.power(m,-1), m], [-k*np.power(m,k), -k*m]])
    [c_1tlambda, c_2tlambda]=np.matmul(M_0i,v_1)
    rhs =r*np.matmul(M_11,[c_1,1-c_2])
    [c_1t, c_2t]=np.matmul(M_0i,rhs)
    w = [c_1*np.power(m,i)+c_2 for i in range(7)]
    e_0 =[c_1t*np.power(m,i)+c_2t +m*r*(1-c_2)*i+r*c_1*i* np.power(m,i)for i in range(7)]
    e_0_lambda =[c_1tlambda*np.power(m,i)+c_2tlambda for i in range(7)]
    rv = nbinom(6,1-p)
    sum_lambda = 1-rv.cdf(5)+np.sum([(e_0_lambda[6-i])*rv.pmf(i-1) for i in range(1,7)])
    sum_no_lambda = np.sum([(e_0[6-i]+i*(w[6-i]-1))*rv.pmf(i-1) for i in range(1,7)])
    ratio = -sum_no_lambda/sum_lambda
    return(ratio)
    
   
   
import pandas as pd
df = pd.DataFrame()
df['hashrate']=[np.float(i)/40 for i in range(1,20)]
df['break-even']=[ratio_six_block_ds(np.float(i)/40) for i in range(1,20)]
df.to_csv("breakeven.csv", index=False)

   
def df_to_latex_table(df):
    f = open("latex_table.tex", "w")
    f.write("\[%")
    f.write("\n")
    f.write('\\begin{tabular}[c]{|')
    col = df.columns
    for c in col: f.write("c|")
    f.write("}\\hline \n")
    for c in col[:-1]: f.write(c + " & ")
    f.write(col[-1])
    f.write( "\\\\\\hline \n")
    for k in df.index[:-1]:
        for a in df.loc[k][:-1]: f.write(str(a)+' & ')
        f.write(str(df.loc[k][-1]))
        f.write( "\\\ \n")
    k = df.index[-1]
    for a in df.loc[k][:-1]: f.write(str(a)+' & ')
    f.write(str(df.loc[k][-1]))
    f.write( "\n")
    f.write("\\\\\\hline \n")
    f.write("\end{tabular} \n")
    f.write("\]")
    f.close()   
        
    
    
    

   
   
### doublespend simulations



def round(l, p):
    state_vector1 = -1
    state_vector2 = 1
    mainchain = 0
    while(mainchain<6):
        if np.random.uniform()<p:
            state_vector2+=1
            state_vector1 += -1 
        else : 
            mainchain+=1
            state_vector1 += 1
    while(1):
        if state_vector1 == -1 : return(l)    
        if state_vector1 == 6 : return(-state_vector2) 
        if np.random.uniform()<p:
            state_vector2+=1
            state_vector1 += -1 
        else : 
            state_vector1 += 1
        
    
    
    
    
    
    
def sim(l,p, N):
    v = 0 
    for i in range(N):
        v = v+round(l,p)
    return(np.float(v)/N)
    
    
def simWL(l,p,N):
    wins = 0
    loss = 0
    for i in range(N): 
        v = round(l,p)
        if v>0:wins+=1
        if v<0: loss+=1
    return(wins,loss)
    
### what win prob should be 

def winprob(p):
    rv = nbinom(6,1-p)
    m = p/(1-p)
    M_0 =np.asarray( [[np.power(m,-1), 1], [np.power(m,k), 1]])
    v_1 = np.asarray([1,0])
    M_0i = np.linalg.inv(M_0)
    [c_1,c_2] = np.matmul(M_0i,v_1)
    w = [c_1*np.power(m,i)+c_2 for i in range(7)]
    total_prob = np.sum([rv.pmf(i)*w[5-i] for i in range(6)])+1-rv.cdf(5)
    return(total_prob)
    








import numpy as np
import scipy.special as sp
from scipy.stats import nbinom


p = 0.3
eps_1 = .01
eps_2 = .01
q = 1- p - eps_1-eps_2
k = 3

lamma = 0.001








#### censorship
def table(lamma, p, eps_1, eps_2, k):
    a = statsNN(p, eps_1, eps_2, k)
    T11 = (a[0]*lamma +a[1], a[2]*lamma + a[3])
    a = statsNC(p, eps_1, eps_2, k)
    T12 = (a[0]*lamma + a[1], a[2]*lamma + a[3])
    a = statsNC(p, eps_2, eps_1, k)
    T21 = (a[2]*lamma + a[3], a[0]*lamma + a[1])
    a = statsCC(p, eps_1, eps_2, k)
    T22 = (a[0]*lamma + a[1], a[2]*lamma + a[3])
    return(T11, T12, T21, T22)
    
   
def normalized_table(lamma, p, eps_1, eps_2, k):
    a = statsNN(p, eps_1, eps_2, k)
    T11 = ((a[0]*lamma +a[1])/eps_1, (a[2]*lamma + a[3])/eps_2)
    a = statsNC(p, eps_1, eps_2, k)
    T12 = ((a[0]*lamma +a[1])/eps_1, (a[2]*lamma + a[3])/eps_2)
    a = statsNC(p, eps_2, eps_1, k)
    T21 = ((a[2]*lamma +a[3])/eps_1, (a[0]*lamma + a[1])/eps_2)
    a = statsCC(p, eps_1, eps_2, k)
    T22 = ((a[0]*lamma +a[1])/eps_1, (a[2]*lamma + a[3])/eps_2)
    return(T11, T12, T21, T22)
    
    
def ntcsv(lamma, p, eps_1, eps_2, k):
    B = normalized_table(lamma, p, eps_1, eps_2, k)
    df = pd.DataFrame(np.asarray(B))
    df.to_csv("table.csv")

import pandas as pd
df = pd.DataFrame(np.asarray(B))

##NN test code
m = p/(1-p)
q = 1- p - eps_1-eps_2 
M_0 =np.asarray( [[np.power(m,-1), 1], [np.power(m,k), 1]])
v_1 = np.asarray([1,0])
M_0i = np.linalg.inv(M_0)
[c_1,c_2] = np.matmul(M_0i,v_1)
M_1 = np.asarray( [[1, -1], [-k*np.power(m,k+1), k]])
w_1 = c_1*m+c_2
unique to pool 1
[ct_1,ct_2] = (eps_1/(1-2*p))*np.matmul(M_0i,M_1).dot([c_1,c_2])
e_1 = ct_1*m+c_2-eps_1*c_2/(1-2*p)+eps_1*c_1*m*m/(1-2*p)
E_1NcNcA = eps_1/(1-p)
E_1NcNcB = -(eps_1/(1-p))*w_1/(1-w_1)-e_1/(1-w_1)
#unique to pool 2
[ct_1,ct_2] = (eps_2/(1-2*p))*np.matmul(M_0i,M_1).dot([c_1,c_2])
e_2 = ct_1*m+c_2-eps_2*c_2/(1-2*p)+eps_2*c_1*m*m/(1-2*p)
E_2NcNcA = eps_2/(1-p)
E_2NcNcB = -(eps_2/(1-p))*w_1/(1-w_1)-e_2/(1-w_1)





def statsNN(p, eps_1, eps_2, k):
    m = p/(1-p)
    q = 1- p - eps_1-eps_2
    M_0 =np.asarray( [[np.power(m,-1), 1], [np.power(m,k), 1]])
    ##NCNC
    v_1 = np.asarray([1,0])
    M_0i = np.linalg.inv(M_0)
    [c_1,c_2] = np.matmul(M_0i,v_1)
    M_1 = np.asarray( [[1, -1], [-k*np.power(m,k+1), k]])
    w_1 = c_1*m+c_2
    #unique to pool 1
    [ct1_1,ct1_2] = (eps_1/(1-2*p))*np.matmul(M_0i,M_1).dot([c_1,c_2])
    e_1 = ct1_1*m+ct1_2-eps_1*c_2/(1-2*p)+eps_1*c_1*m*m/(1-2*p)
    E_1NcNcA = eps_1/(1-p)
    E_1NcNcB = -(eps_1/(1-p))*w_1/(1-w_1)-e_1/(1-w_1)
    #unique to pool 2
    [ct2_1,ct2_2] = (eps_2/(1-2*p))*np.matmul(M_0i,M_1).dot([c_1,c_2])
    e_2 = ct2_1*m+ct2_2-eps_2*c_2/(1-2*p)+eps_2*c_1*m*m/(1-2*p)
    E_2NcNcA = eps_2/(1-p)
    E_2NcNcB = -(eps_2/(1-p))*w_1/(1-w_1)-e_2/(1-w_1)
    return(E_1NcNcA,E_1NcNcB,E_2NcNcA,E_2NcNcB)

def statsNC(p, eps_1, eps_2, k):
    m = p/(1-p)
    q = 1- p - eps_1-eps_2
    pt = p+eps_2
    M3 = np.asarray( [[1, 1,-1], [m, 1, -1/(1-pt)],[np.power(m,k), 1,0]])
    v_3 =  np.asarray([0,-pt/(1-pt),0])
    M3i = np.linalg.inv(M3)
    [c_1,c_2, b] = np.matmul(M3i,v_3)
    w_1 = c_1*m+c_2
    #Pool1
    M4 =  np.asarray( [[0, 0,0], [-w_1/(1-pt), -m*m/(1-2*p), 1/(1-2*p)],[0,-k*np.power(m,k+1)/(1-2*p), k/(1-2*p)]])
    [ct1_1,ct1_2, b_1] = eps_1*M3i.dot(M4).dot([1,c_1,c_2])
    e_1 = (b_1 - eps_1*w_1)/(1-pt)
    E_1NcCA = eps_1/(1-pt)
    E_1NcCB = -(eps_1/(1-pt))*w_1/(1-w_1)-e_1/(1-w_1)
    #Pool2
    M4b =  np.asarray( [[0, 0], [ -m*m/(1-2*p), 1/(1-2*p)],[-k*np.power(m,k+1)/(1-2*p), k/(1-2*p)]])
    [ct2_1,ct2_2, b_2] = eps_2*M3i.dot(M4b).dot([c_1,c_2])
    e_2 = b_2/(1-pt)
    E_2NcCA = 0
    E_2NcCB = -e_2/(1-w_1)
    return(E_1NcCA,E_1NcCB,E_2NcCA,E_2NcCB)

def statsCC(p, eps_1, eps_2, k):
    m = p/(1-p)
    q = 1- p - eps_1-eps_2
    pt = p+eps_1+eps_2
    M3 = np.asarray( [[1, 1,-1], [m, 1, -1/(1-pt)],[np.power(m,k), 1,0]])
    v_3 =  np.asarray([0,-pt/(1-pt),0])
    M3i = np.linalg.inv(M3)
    [c_1,c_2, b] = np.matmul(M3i,v_3)
    M4b =  np.asarray( [[0, 0], [ -m*m/(1-2*p), 1/(1-2*p)],[-k*np.power(m,k+1)/(1-2*p), k/(1-2*p)]])
    w_1 = c_1*m+c_2
    #Poo1 
    [ct1_1,ct1_2, bt_1] = eps_1*M3i.dot(M4b).dot([c_1,c_2])
    e_1 = bt_1/(1-pt)
    E_1CCA = 0
    E_1CCB =-e_1/(1-w_1)
    #Poo1 2
    [ct2_1,ct2_2, bt_2] = eps_2*M3i.dot(M4b).dot([c_1,c_2])
    e_2 = bt_2/(1-pt)
    E_2CCA = 0
    E_2CCB =-e_2/(1-w_1)
    return(E_1CCA,E_1CCB,E_2CCA,E_2CCB)





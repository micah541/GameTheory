

from itertools import permutations
import numpy as np


def get_cut(q,v):
    sum = 0.0
    for j in q:
        sum=sum+v[j]
        if sum>0.5: return(j)

def Shapley(w):
    k = len(w)
    perm = permutations(range(k))
    wins = np.zeros(k)
    for p in perm:wins[get_cut(p,w)]+=1
    wins = wins/np.math.factorial(k)
    return(wins)


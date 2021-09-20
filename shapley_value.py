from itertools import permutations

w = [.25, .25, .3, .1, .1]

k = len(w)


for i in perm:
    sum = 0
    for j in range(k):
        sum = sum+i[j]
        if sum>0.5 : return()

def get_cut(p, w):
    sum=0
    for j in range(k):
        sum = sum+w[i[j]]
        if sum>0.5 : return(j)

perm = permutations(range(k))
wins = np.zeros(k)
for p in perm:
    wins[get_cut(p,w)]+=1


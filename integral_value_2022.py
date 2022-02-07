#The file titled "integral_value_2022.py" contains the python code for the graphical algorithm for the integration of monomials in the Chow ring of the moduli space of stable marked curves of genus zero. The algorithm is called "the forest algorithm" in the preprint arXiv:2102.03575. 

import itertools
from collections import Counter
from ast import literal_eval
import networkx as nx
from itertools import chain, combinations

def monomial_value(m): # m is a monomial in the ambient Chow ring; output the integral value of this monomial. See Section 1 of the preprint arXiv:2102.03575 for related background and the definition of the integral value of a monomial in the ambient Chow ring.
    z = set(m[0][0]).union(set(m[0][1]))
    if len(m) != len(z)-3:
        return 0
    else:
        if keel(m) == 0: 
            return 0
        else:
            t = redundancy_tree(m)
            s = sign(m)
            a = absolute_value(t)
            return s*a
  

def keel(m):  # m is a monomial in the ambient Chow ring; output 1 if m is a tree monomial, output 0 otherwise. See Section 2 of the preprint arXiv:2102.03575 for the definition of a tree monomial.
    for i in range(len(m)):
        for j in range(len(m)):
            if keel_micro(m[i], m[j]) == 0:
                return 0
                break
    return 1 


def keel_micro(m1,m2): #m1,m2 are two tuples representing two elements in the ambient Chow group of codimension one, respectively. Output 0 if these two monomials fulfill the Keel's quadratic relation, output 1 otherwise. See Section 2 of the preprint arXiv:2102.03575 for the definition of Keel's quadratic relation.
    if (set(m1[0]).intersection(set(m2[0])) != set() and set(m1[0]).intersection(set(m2[1])) != set() and set(m1[1]).intersection(set(m2[0])) != set() and set(m1[1]).intersection(set(m2[1])) != set()):
        return 0
    else:
        return 1


           


def sign(m):# m is a monomial in the ambient Chow ring; output the sign of the integral value of this monomial. See Section 3 of the preprint arXiv:2102.03575 for the calculation of the sign of the integral value of the input monomial.
    a = 0
    e = Counter(m)
    l = len(e)
    for x in e:
        a = a + e[x]
    return pow(-1,a-l)



def redundancy_tree(m):# m is a monomial in the ambient Chow ring; output the redundancy tree of this monomial. See Section 3 of the preprint arXiv:2102.03575 for the definition of the redundancy tree of the input monomial.
    G = nx.DiGraph()
    G1 = nx.Graph()
    P = []
    C = []
    x = m[0][0]
    y = m[0][1]
    z = set(x).union(set(y))
    w = m[0]
    del m[0]
    for i in m:
        P.append(i[0])
        P.append(i[1])
    for j in P:
        if (set(j).issubset(set(x)) or set(j).issubset(set(y))):
            C.append(j)
    C.append(x)
    C.append(y)
    G.add_edge(x,y)
    for c in C:
        G.add_node(c)
        for d in C:
            if ( d != c and set(c).issubset(set(d))):
                G.add_edge(c, d)
    L = list(G.edges)   
    for e in L:
        for f in L:
            if ( f != e and f[0] == e[0]):
                if set(f[1]).issubset(set(e[1])):
                    L.remove(e)
                else:
                    if set(e[1]).issubset(set(f[1])):
                        L.remove(f)
    G1.add_nodes_from(C)                    
    G1.add_edges_from(L)
    for c in C:
        c1 = c
        for n in G1[c]:
            if set(n).issubset(set(c)):
                c = tuple(set(c).difference(set(n)))
        a = len(c) + G1.degree(c1) - 3
        G1.node[c1]['weight'] = a
    G1.remove_edges_from(L)
    m.insert(0,w)
    for l in L:
        G1.add_node(l)
        l1 = tuple(z.difference(set(l[0])))
        w1 = tuple((l[0],l1))
        w2 = tuple((l1, l[0]))
        if m.count(w1) != 0:
            b = m.count(w1)
        else:
            b = m.count(w2)   
        G1.node[l]['weight'] = b-1
        G1.add_edge(l,l[0])
        G1.add_edge(l,l[1])
    return G1


def binomialCoeff(n, k):#This function outputs the binomial coefficient of n choose k.
    if n<k:
        result=0
    else:
        result = 1
        for i in range(1, k+1):
            result = result * (n-i+1) / i
    return result        
  
  
def absolute_value(g):#g is a redundancy tree. Input is a redundancy tree, output is the value of the redundancy tree, i.e. the absolute value of the monomial m if g is the redundancy tree of monomial m. See Section 3 of the preprint arXiv:2102.03575 for the definition of the redundancy tree of the input monomial and that of the value of a redundancy tree.
    g1 = g.nodes   
    g2 = [x for x in g.nodes() if g.degree(x)==0]
    g3 = [x for x in g.nodes() if g.degree(x)==1]
    if g2 != []:
        return 0
    elif g3 == []:
        return 1
    else:
        for x in g3:
            s0 = g1[x]['weight']                
            s1 = list(g.neighbors(x))[0]
            s2 = g1[s1]['weight']
            if s0 == s2:
                #print 121
                g.remove_node(x)
                g.remove_node(s1)
                var = absolute_value(g)
                return var
            elif s0 > s2:
                return 0
            else:    
                #print 123
                g1[s1]['weight'] = s2-s0    
                g.remove_node(x)                
                var = binomialCoeff(s2,s0)*absolute_value(g)        
                return var  
  
  
#In the sequel are three examples showing what form should a monomial be as the input of the function "monomial_value".  
  
m = [ ((2,5),(1,3,4,6,7,8)), ((2,5),(1,3,4,6,7,8)),((3,4,8),(1,2,5,6,7)), ((3,8),(1,2,4,5,6,7)), ((2,5,6),(1,3,4,7,8)), ((1,7),(2,3,4,5,6,8))]
#This is the encoding of the monomial \delta^2_{25|134678}\cdot \delta_{348|12567}\cdot \delta_{38|124567}\cdot \delta_{256|13478}\cdot \delta_{17|234568}.

m1 = [((1,2),(3,4,5,6)), ((1,2),(3,4,5,6)), ((1,2,3,4),(5,6))]
#This is the encoding of the monomial \delta^2_{12|3456}\cdot \delta_{1234|56}.


m2 = [((1,2,3),(4,5,6,7,8,9,10,11,12,13,14)), ((1,2,3),(4,5,6,7,8,9,10,11,12,13,14)), ((1,2,3),(4,5,6,7,8,9,10,11,12,13,14)), ((9,10,11),(1,2,3,4,5,6,7,8,12,13,14)), ((9,10,11),(1,2,3,4,5,6,7,8,12,13,14)),((9,10,11),(1,2,3,4,5,6,7,8,12,13,14)),((9,10,11),(1,2,3,4,5,6,7,8,12,13,14)),((9,10,11),(1,2,3,4,5,6,7,8,12,13,14)), ((8,12,13,14), (1,2,3,4,5,6,7,9,10,11)), ((12,13,14),(1,2,3,4,5,6,7,8,9,10,11)), ((12,13,14),(1,2,3,4,5,6,7,8,9,10,11))]

print monomial_value(m2)# Here is one example of how to apply the function "monomial_value" to the monomial m2 given above. This coincides with Example 3.1 of the preprint arXiv:2102.03575.


from solver import *
from vecutil import *
from matutil import *

def rep2vec(u, veclist):
    M = rowdict2mat(veclist)
    return solve(M,u)

def vec2rep(veclist, v):
    M = coldict2mat(veclist)
    return solve(M,v)

def is_superfluous(L, i):
    if len(L) <= 1 : return False
    A = coldict2mat(L[0:i]+L[i+1:len(L)])
    b = L[i]
    u = solve(A,b)
    residual = b - A*u
    return True if residual * residual < 1.0/(10**14) else False

def is_independent(L):
    r = [is_superfluous(L,i) for i in range(len(L))]
    return True if sum(r) == 0 else False

def subset_basis(T):
    #if len(T) == 1: return T
    S = []
    for v in T:
        if is_independent(S+[v]): S = S + [v]
    return S

def superset_basis(T, L):
    S = T
    for v in L:
        if is_independent(S+[v]): S = S + [v]
    return S

def exchange(S, A, z):
    B = [v for v in S if v not in A]
    beta = vec2rep(A+B, z)
    #beta = solve(coldict2mat(A+B), z)
    Alen = len(A)
    for i in range(len(B)):
        j = Alen + i
        if beta[j] != 0:
            betaj = beta[j]
            beta[j] = 0
            return 1/betaj*z -1/betaj*coldict2mat(A+B)*beta
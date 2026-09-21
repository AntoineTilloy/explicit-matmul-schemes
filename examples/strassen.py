"""Run the package's exact primitive reader on a manageable example."""
from fractions import Fraction
from cwschemes._resources import index, materialize
for name in index()['paths']:
    if name.startswith('seeds/'):materialize(name)
from cwschemes._engine.general_square import Ordinary

s=Ordinary({'kind':'lille','n':2,'rank':7})
A=[[1,2],[3,4]];B=[[5,6],[7,8]]
C=[[Fraction(0) for _ in range(2)] for _ in range(2)]
for t in range(s.rank):
    a=sum(s.factor_entry(0,t,i,j)*A[i][j] for i in range(2) for j in range(2))
    b=sum(s.factor_entry(1,t,i,j)*B[i][j] for i in range(2) for j in range(2))
    for i in range(2):
        for j in range(2):C[i][j]+=s.factor_entry(2,t,i,j)*a*b
assert C==[[19,22],[43,50]]
print(C)

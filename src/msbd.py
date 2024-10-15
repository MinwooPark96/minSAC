from .model import CausalDiagram as CD
from typing import Iterable, Optional, Union, FrozenSet
from src.order import *
from src.basis import SequentialAdjustment

class MultiOutcomeSequentialBackdoorCriterion(SequentialAdjustment):

    def __init__(self, 
                 G: CD, 
                 xs: Iterable[str],
                 ys: Iterable[str],
                 zs: Iterable[str] = frozenset()):
        
        self.G = G
        self.V = G.V
        self.xs, self.ys, self.zs = wrap(xs), wrap(ys), wrap(zs)
        
        assert self.xs & self.ys & self.zs == set()
        
        self.X = OrderedX(self.G, self.xs)
        self.Y = OrderedY(self.G, self.X, ys)
        self.Z = Orderedmsbd(self.G,self.X,zs)
        self.m = self.X.m
        
    def localGraph(self,i: int) -> CD:
        return self.G.underbar(self.X[i]).do(self.X.geq(i+1))    
    
    def condition1(self,i: int) -> bool:
        return all(zi not in self.G.De(self.X.geq(i)) for zi in self.Z[i])
    
    def condition2(self,i: int) -> bool:
        L = self.localGraph(i)
        return L.independent(self.X[i],self.Y.geq(i),self.H(i-1)|self.Z[i]|self.X.geq(i+1))
        
    def criteria(self) -> bool:
        for i in range(1,self.m+1):
            if not self.condition1(i):
                print(f"condition1 (descendant) is violated at i={i}")
                return False
            if not self.condition2(i):
                print(f"condition2 (independent) is violated at i={i}")
                return False
        else:
            return True 
    
    def __getitem__(self,i:int):
        assert i>= 1 and i <= self.m
        return self.localGraph(i)
    
    def __repr__(self):
        return f'MultiOutcomeSequentialBackdoorCriterion({self.G},\n{self.X},\n{self.Y},\n{self.Z})'
    
    def __str__(self):
        return f'MultiOutcomeSequentialBackdoorCriterion({self.G},\n{self.X},\n{self.Y},\n{self.Z})'

mSBD = MultiOutcomeSequentialBackdoorCriterion
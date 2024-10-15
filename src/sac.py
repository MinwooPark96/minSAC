from .model import CausalDiagram as CD
from .utils import fzset_union, optional_next
from typing import Iterable, Optional
from src.order import *
from src.basis import SequentialAdjustment
from src.ac import AdjustmentCriterion
from . import navigator

class SequentialAdjustmentCriterion(SequentialAdjustment):
    
    def __init__(self, 
                 G: CD, 
                 xs: Iterable[str], 
                 ys: Iterable[str], 
                 zs: Iterable[str] = frozenset()) -> None: 
        
        self.G = G
        self.V = G.V
        self.xs, self.ys, self.zs = wrap(xs), wrap(ys), wrap(zs)
        assert self.xs & self.ys & self.zs == set()
        
        self.X = OrderedX(self.G, self.xs) 
        self.order_iter = iter(self.X.get_all_causal_order())
        optional_next(self.order_iter)
        
        self.Y = OrderedY(self.G,self.X,ys)
        self.Z = Orderedsac(self.G,self.X,zs)
        self.m = self.X.m 
        
    def __next(self) -> Optional[list[str]]:
        next_order = optional_next(self.order_iter)
        return next_order

    def set_next(self) -> bool:
        next_order = self.__next()
        if next_order:
            self.X.set_partition(next_order)
            return True
        else :
            return False
    
    def dpcpSet(self, i: int):
        W = navigator.dpcpSet(self.G,self.X[i],self.Y.geq(i))
        return self.G.De(W)
    
    def condition1(self, i:int) -> bool:
        return not (self.Z[i] & self.dpcpSet(i))
    
    def condition2(self, i:int) -> bool:    
        return bool(self.psbd(i).independent(self.X[i],self.Y.geq(i), self.H(i-1) | self.Z[i]))
    
    def criteria(self):
        
        for i in range(1,self.m+1):
            if not self.condition1(i):
                print(f"condition1 (z | dpcp) is violated at i={i}")
                return False
            if not self.condition2(i):
                print(f"condition2 (opend psbd) is violated at i={i}")
                return False
        else:
            return True 
    
    def psbd(self,i:int):
        return navigator.pbd(self.G.do(self.X.geq(i+1)),self.X[i],self.Y.geq(i))
    
    def F(self, i: int):
        return self.X | self.Y | self.H(i-1) | self.dpcpSet(i) | self.G.De(self.X.geq(i+1))
    
    def construct(self):
        for i in range(1,self.m+1): # 1,2,...,m
            psbd_i = self.psbd(i)
            zi = psbd_i.An(self.X[i]|self.Y.geq(i)|self.H(i-1)) - self.F(i)
            self.set_covariate(i,zi)
        
    def minsac(self):
        self.construct() #self.Z <- Z(an)
        minz = [frozenset()]
        
        for i in range(1,self.m+1): # 1,2,...,m
            Zian = self.Z[i]
            Zi1 = navigator.closure(self.G,self.Y.geq(i),self.X[i],Zian | self.H(i-1)) & Zian
            Zimin = navigator.closure(self.G,self.X[i],self.Y.geq(i),Zi1 | self.H(i-1)) & Zi1
            minz.append(wrap(Zimin))
            
        for i in range(1,self.m+1):
            self.set_covariate(i,minz[i])
        
    def to_ac(self) -> AdjustmentCriterion:
        subZ = fzset_union([self.Z[i] for i in range(1,self.m+1)])
        return AdjustmentCriterion(G = self.G,
                                   xs = self.xs,
                                   ys = self.ys,
                                   zs = subZ)
    
    def __getitem__(self,i:int):
        assert i>= 1 and i <= self.m
        return self.psbd(i)    

    def __repr__(self):
        return f'SequentialAdjustmentCriterion({self.G},\n{self.X},\n{self.Y},\n{self.Z})'
    
    def __str__(self):
        return f'SequentialAdjustmentCriterion({self.G},\n{self.X},\n{self.Y},\n{self.Z})'

    
SAC = SequentialAdjustmentCriterion


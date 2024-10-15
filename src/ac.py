from .model import CausalDiagram as CD
from .utils import wrap
from . import navigator
from typing import Iterable

import itertools


class AdjustmentCriterion:
    
    def __init__(self, 
                 G: CD, 
                 xs: Iterable[str], 
                 ys: Iterable[str], 
                 zs :Iterable[str] = frozenset()):
        
        self.G = G
        self.V = G.V        
        self.xs, self.ys = wrap(xs), wrap(ys)
        self.F = navigator.dpcpSet(G = self.G,
                                   xs = self.xs,
                                   ys = self.ys)
        self.z0 = self.explicit_admissible_set()
        self.zs = frozenset()
        self.set_covariate(zs)
        
    def set_covariate(self, zs: Iterable[str]):
        zs = wrap(zs)
        self.zs = zs
    
    def get_all_covariates(self, with_empty_set: bool = False) -> list[frozenset[str]]: 
        Z = self.G.V - self.xs - self.ys - self.F 
        G_pbd = navigator.pbd(self.G,self.xs,self.ys)
        
        covariates = list()
        
        min_size = 0 if with_empty_set else 1
        
        for size in range(min_size,len(Z)+1): # size 1 로 수정
            for subset in itertools.combinations(Z, size):
                if G_pbd.independent(
                    xs = self.xs,
                    ys = self.ys,
                    zs = subset):
                    covariates.append(wrap(subset))
        return covariates
    
    def get_all_backdoor_covariates(self, with_empty_set: bool = False) -> list[frozenset[str]]: 
        Z = self.G.V - self.G.De(self.xs) - self.ys
        covariates = list()
        min_size = 0 if with_empty_set else 1

        for size in range(min_size,len(Z)+1): # size 1 로 수정
            for subset in itertools.combinations(Z, size):
                if self.G.underbar(self.xs).independent(
                    xs = self.xs,
                    ys = self.ys,
                    zs = subset):            
                    covariates.append(wrap(subset))
        return covariates
    
    def explicit_admissible_set(self) -> frozenset[str] :
        z0 = self.G.An(self.xs|self.ys)- self.xs - self.ys - self.F
        return z0
    
    def criteria(self) :
        return self.zs in self.get_all_covariates(with_empty_set = True)
        
    def __repr__(self):
        return f'Adjustment Criterion({self.G},\n{self.xs},\n{self.ys}),\n{self.zs})'
    
    def __str__(self):
        return f'Adjustment Criterion({self.G},\n{self.xs},\n{self.ys}),\n{self.zs})'
    


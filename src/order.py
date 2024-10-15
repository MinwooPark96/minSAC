from .model import CausalDiagram as CD
from .utils import fzset_union, wrap
from typing import Iterable, AbstractSet, Union
from abc import ABC
import networkx as nx
import functools

class Orderedset(AbstractSet, ABC):
    
    def leq(self, i: int) -> frozenset[str]:
        if i <= 0:
            return frozenset()
        elif i > self.m:
            return fzset_union([self[j] for j in range(1,self.m+1)])
        return fzset_union([self[j] for j in range(1,i+1)])
        
    
    def geq(self, i: int) -> frozenset[str]:
        if i <= 0:
            return fzset_union([self[j] for j in range(1,self.m+1)])
        elif i > self.m:
            return frozenset()
        return fzset_union([self[j] for j in range(i,self.m+1)])
        
    def __contains__(self, x: object) -> bool:
        return x in self.V
    
    def __iter__(self):
        return iter(self.V)
    
    def __next__(self):
        return next(self.V)
    
    def __len__(self):
        return self.m

    # __getitem__ shold be overrided in OrderedY 
    def __getitem__(self,i):        
        return  self.P[i-1] if 0 < i <= self.m else frozenset()
        
    def __str__(self):
        return str(self.P)
    
    def __repr__(self):
        return str(self.P)
    
    def __or__(self, other):
        # if isinstance(other,OrderedX) or isinstance(other,OrderedY) or isinstance(other,Orderedmsbd):
        if hasattr(other,'V'):
            return self.V | other.V
        else :
            return self.V | wrap(other)
    
    def __and__(self, other):
        # if isinstance(other,Orderedset):
        if hasattr(other,'V'):
            return self.V | other.V
        else :
            return self.V | wrap(other)
    
    def __eq__(self, other):
        return self.V == other.V
    
class OrderedX(Orderedset):        
    
    def __init__(self,
                 G: CD,
                 xs: Iterable[str]):
        
        self.G = G
        self.V= wrap(xs)
        self.get_all_causal_order = functools.lru_cache()(self.get_all_causal_order)
        self.P= self.partition()
        self.m = len(self.V)
        
    def get_all_causal_order(self) -> list[list[str]]: 
        LP = self.G.LP(self.V)
        gg = nx.DiGraph(LP.edges)
        gg.add_nodes_from(self.V)
        top_to_bottom = list(nx.all_topological_sorts(gg))
        return top_to_bottom
        
    def partition(self) -> list[frozenset[str]]:
        result = list()
        for x in self.get_all_causal_order()[0]:
            if x in self.V:
                result.append(wrap(x))
        return result
    
    def set_partition(self, xs : Iterable[str]):
        xs = wrap(xs)
        assert xs == self.V
        for i in range(0,self.m):
            self.P[i] = wrap(xs[i])

class OrderedY(Orderedset):        
    def __init__(self,
                 G: CD,
                 X: OrderedX,
                 ys: Iterable[str] = frozenset()):
        self.G = G
        self.X = X
        self.V = self.ys = wrap(ys)
        self.m = len(self.X)
        self.P = self.Y = self.partition()
            
    def partition(self) -> list[frozenset[str]]:
        Y0 = self.ys - self.G.De(self.X.V)
        ps = [Y0]
        for i in range(1,self.m+1):
            left = self.V - fzset_union(ps)
            right = self.G.De(self.X[i]) - self.G.De(self.X.geq(i+1))
            ps.append(left & right)
        return ps
    
    def __getitem__(self,i):        
        return self.P[i] if 0 <= i <= self.m else frozenset()
    
    def leq(self,i:int):
        if i < 0:
            return frozenset()
        elif i > self.m:
            return fzset_union([self[j] for j in range(0,self.m+1)])
        return fzset_union([self[j] for j in range(0,i+1)])
        
    def geq(self,i:int):
        if i < 0:
            return fzset_union([self[j] for j in range(0,self.m+1)])
        elif i > self.m:
            return frozenset()
        return fzset_union([self[j] for j in range(i,self.m+1)])
    
    
class Orderedmsbd(Orderedset):
    def __init__(self,
                 G: CD,
                 X: OrderedX,
                 zs: Iterable[str] = frozenset()):
        
        self.G = G
        self.X = X 
        self.V = wrap(zs)
        self.m = self.X.m
        self.P = self.partition()
    
    def set_partition(self, i: int, zi: Iterable[str]) :
        zi = wrap(zi)
        assert i > 0 and i <= self.m
        self.P[i-1] = zi
    
    def partition(self) -> list[frozenset[str]]:
        ps = list()
        Z = self.V    
        for i in range(1,self.m+1):
            Zi = list()
            for v in self.V:
                if v in Z - self.G.de(self.X.geq(i+1)): #search space
                    Zi.append(v)          
            if Zi:
                ps.append(wrap(Zi))
            else :
                ps.append(frozenset())
            Z = Z - wrap(Zi)
        return ps

class Orderedsac(Orderedset):   

    def __init__(self,
                 G: CD,
                 X: OrderedX,
                 zs: Iterable[str] = frozenset()
                 ):
        
        self.G = G
        self.X = X 
        self.V = wrap(zs)
        self.m = self.X.m
        self.P = self.partition()
        
    def set_partition(self,i: int, zi: Iterable[str]):
        zi = wrap(zi)
        assert i > 0 and i <= self.m
        assert all(z in self.V for z in zi)
        self.P[i-1] = wrap(zi)
    
    def partition(self) -> list[frozenset[str]]:
        pset = list()
        Z = self.V
        for i in range(1,self.m+1):
            zi = list()
            for v in self.V:
                if v in (Z - self.G.De(self.X.geq(i))): #search space
                    zi.append(v)         
            if zi:
                pset.append(wrap(zi))
            else :
                pset.append(frozenset())
            
            Z = Z - wrap(zi)
        
        return pset
        

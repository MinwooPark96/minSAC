from .model import CausalDiagram as CD
from .utils import fzset_union, wrap
from typing import Iterable, AbstractSet, Union
from abc import ABC
import networkx as nx

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
            return self.V & other.V
        else :
            return self.V & wrap(other)
    
    def __eq__(self, other):
        return self.V == other.V
    
class OrderedX(Orderedset):        
    
    def __init__(self,
                 G: CD,
                 xs: Iterable[str]):
        
        self.G = G
        self.V = wrap(xs)
        if not self.V or not self.V <= G.V:
            raise ValueError("Treatments must be a nonempty subset of graph vertices")
        directed = nx.DiGraph()
        directed.add_nodes_from(sorted(G.V))
        directed.add_edges_from(G.edges)
        if not nx.is_directed_acyclic_graph(directed):
            raise ValueError("The directed part of the causal graph must be acyclic")
        self.m = len(self.V)
        # One deterministic order; do not enumerate all orders during construction.
        order = [v for v in nx.lexicographical_topological_sort(directed) if v in self.V]
        self.P = [wrap(v) for v in order]
        
    def get_all_causal_order(self) -> list[list[str]]: 
        """Explicitly enumerate all orders (potentially factorial output)."""
        return list(self.iter_causal_orders())

    def iter_causal_orders(self):
        """Lazily enumerate orders respecting paths through other vertices."""
        gg = nx.DiGraph()
        gg.add_nodes_from(sorted(self.V))
        for x in sorted(self.V):
            gg.add_edges_from((x, y) for y in sorted((self.G.De(x) & self.V) - {x}))
        return nx.all_topological_sorts(gg)
        
    def partition(self) -> list[frozenset[str]]:
        return list(self.P)
    
    def set_partition(self, xs : Iterable[str]):
        xs = [xs] if isinstance(xs, str) else list(xs)
        if len(xs) != self.m or frozenset(xs) != self.V:
            raise ValueError("Treatment order must contain each treatment exactly once")
        for i, x in enumerate(xs):
            if self.G.De(x) & set(xs[:i]):
                raise ValueError("Treatment order must respect directed ancestry")
        self.P = [wrap(x) for x in xs]

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
        if not 1 <= i <= self.m:
            raise ValueError("Covariate stage is outside the treatment order")
        if not zi <= self.G.V - self.X.V:
            raise ValueError("Covariates must be observed nontreatment vertices")
        if zi & self.G.De(self.X.geq(i + 1)):
            raise ValueError("Covariates cannot descend from a later treatment")
        if any(zi & block for j, block in enumerate(self.P, 1) if j != i):
            raise ValueError("Covariate blocks must be pairwise disjoint")
        self.P[i-1] = zi
        self.V = fzset_union(self.P)
    
    def partition(self) -> list[frozenset[str]]:
        pset = list()
        Z = self.V
        for i in range(1,self.m+1):
            zi = list()
            for v in self.V:
                if v in (Z - self.G.De(self.X.geq(i+1))): # current descendants are allowed
                    zi.append(v)         
            if zi:
                pset.append(wrap(zi))
            else :
                pset.append(frozenset())
            
            Z = Z - wrap(zi)
        
        return pset
        

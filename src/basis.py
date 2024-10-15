from .utils import wrap
from typing import Iterable
from abc import ABC

class SequentialAdjustment(ABC):
    
    def set_covariate(self, i: int, zi: Iterable[str] = frozenset()):
        zi = wrap(zi)
        self.Z.set_partition(i,zi)
        
    def H(self, i: int) -> frozenset[str]:
        return self.X.leq(i) | self.Y.leq(i) | self.Z.leq(i)
    
    def causalEffect(self):
        raise NotImplementedError
            
    def formula(self): 
        raise NotImplementedError
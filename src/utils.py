from typing import Iterable, Optional, Union, Mapping, Tuple, TypeVar, Generator, Collection, Iterator
from collections import defaultdict
import itertools
from itertools import chain
import networkx as nx
import numpy as np
import sys

VSet = frozenset[str]
VVs = Union[str, Iterable[str]]

T = TypeVar('T')
KT = TypeVar('KT')
KT2 = TypeVar('KT2')
VT = TypeVar('VT')

def shuffled(xs: Iterable[T], reproducible=False) -> list[T]:
    """ A list with shuffled elements """
    if reproducible:
        xs = sorted(xs)
    else:
        xs = list(xs)
    np.random.shuffle(xs)
    return xs


def set_union(sets: Iterable[Collection[T]]) -> set[T]:
    return set(chain(*sets))

def optional_next(xs: Iterator[T], default: Optional[T] = None) -> T:
    try:
        x = next(xs)
        return x
    except StopIteration:
        return default


def default_P_U(mu: Mapping):

    def P_U(d):
        p_val = 1.0
        for k in mu.keys():
            p_val *= (1 - mu[k]) if d[k] == 0 else mu[k]
        return p_val

    return P_U


def as_settuple(xxs:Iterable[Iterable[T]]) -> set[tuple[T]]:
    result = set()
    for xs in xxs:
        result.add(tuple(xs))
    return result

def as_setfzset(xxs:Iterable[Iterable[T]]) -> set[frozenset[T]]:
    result = set()
    for xs in xxs:
        result.add(frozenset(xs))
    return result

def as_tuplefzset(xxs:Iterable[Iterable[T]]) -> tuple[frozenset[T]]:
    return tuple(as_setfzset(xxs))


def as_sortups(xss: Iterable[Iterable[T]]) -> Generator[Tuple[T, ...], None, None]:
    """ Generate tuple-applied elements of a given iterable """
    for xs in xss:
        yield sortup(xs)

def sortup(xs: Iterable[T]) -> Tuple[T, ...]:
    """ Syntactic sugar for tuple(sorted(...)) """
    # sorted tuple
    return tuple(sorted(xs))

def sortup2(xxs: Iterable[Iterable[T]]) -> Tuple[Tuple[T, ...], ...]:
    # twice sorted tuples
    return sortup([sortup(xs) for xs in xxs])

def dict_only(a_dict: dict[KT, VT], keys_only_in: Collection[KT]) -> dict[KT, VT]:
    """ Copy of dictionary with the intersection of the original keys and the specified keys """
    return {k: a_dict[k] for k in keys_only_in if k in a_dict}


def dict_except(a_dict: dict[KT, VT], keys_except: Collection[KT]) -> dict[KT, VT]:
    return {k: v for k, v in a_dict.items() if k not in keys_except}

def fzset_union(sets: Iterable[Collection[T]]) -> frozenset[T]:
    return frozenset(itertools.chain(*sets))

def unique(vs: Iterable[T]) -> list[T]:
    seen = set()
    at = list()
    for v in vs:
        if v not in seen:
            at.append(v)
            seen.add(v)
    return at

def combinations(xs: Iterable[T], reverse = False) -> Generator[Tuple[T, ...], None, None]:
    xs = list(xs)
    if reverse:
        for i in range(len(xs), -1, -1):  # 순서를 반대로 바꿈
            for comb in itertools.combinations(xs, i):
                yield comb
    else:
        for i in range(len(xs) + 1):
            for comb in itertools.combinations(xs, i):
                yield comb

def sortup(xs: Iterable[T]) -> Tuple[T, ...]:
    return tuple(sorted(xs))

def sortup2(xxs: Iterable[Iterable[T]]) -> Tuple[Tuple[T, ...], ...]:
    return sortup([sortup(xs) for xs in xxs])

def remove_reversed_pairs(pairs) -> list[frozenset[str,str]]:
    seen = set()
    result = []
    
    for pair in pairs:
        sorted_pair = sortup(pair)
        if sorted_pair not in seen:
            seen.add(sorted_pair)
            result.append(pair)
    
    result = [frozenset(pair) for pair in result]
    return result

def count_keys(d: Mapping[str, VSet]) -> int:
    """[minwoo] count the number of keys in the dictionary"""
    counter = 0
    for k in d:
        counter += len(d[k])
    return counter

def pairs2dict(xys, backward=False):
    """
    [minwoo] 
        original version is in npsem.model. 
        It is useful to store parent and child information
    """
    dd = defaultdict(set)
    if backward:
        for x, y in xys:
            dd[y].add(x)
    else:
        for x, y in xys:
            dd[x].add(y)

    return defaultdict(frozenset, {key: frozenset(vals) for key, vals in dd.items()})

def pairs2eachdict(xys, is_triple = False):
    dd = defaultdict(set)
    if is_triple: # mainly used for bidirected edges
        for x, y, _ in xys:
            dd[x].add(y)
            dd[y].add(x)
    else :
        for x, y in xys:
            dd[x].add(y)
            dd[y].add(x)
    
    return defaultdict(frozenset, {key: frozenset(vals) for key, vals in dd.items()})

def wrap(v_or_vs: Union[str, Iterable[str]]) -> Optional[VSet]:
    if v_or_vs is None:
        return None
    if isinstance(v_or_vs, str):
        return frozenset({v_or_vs})
    else:
        return frozenset(v_or_vs)

def get_current_function():
    return sys._getframe(1).f_code.co_name

def print_helper(vs:VVs):
        if not vs:
            return ' '
        vs = sorted([v.lower() for v in vs])
        s = ''
        for i in range(len(vs)-1):
            if len(vs[i]) != 1:
                s = s + vs[i][0]+'_'+vs[i][1] + ','
            else :
                s = s + vs[i] + ','
        if len(vs[-1]) != 1:
            s = s + vs[-1][0]+'_'+vs[-1][1]
        else :
            s = s + vs[-1]
        return s

def all_subsets(s):
    """Return a list of all subsets of a set s."""
    from itertools import chain, combinations
    return list(chain.from_iterable(combinations(s, r) for r in range(len(s)+1)))

def subsets_with_complements(s):
    """Generate all subsets of set s and their complements, including the empty set."""
    s = set(s)  # Ensure the input is a set
    subsets = all_subsets(s)
    divisions = []
    
    for subset in subsets:
        subset1 = set(subset)
        subset2 = s - subset1
        divisions.append((subset1, subset2))
    
    return divisions

def with_default(x, dflt=None):
    return x if x is not None else dflt


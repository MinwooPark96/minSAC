from .model import CausalDiagram as CD
from .utils import wrap, set_union
from copy import deepcopy
import itertools
from typing import Iterable

def connected_path(G: CD, x: str, y: str) -> list[list[str]]:
    
    assert x in G.V and y in G.V
    
    if x == y:
        return [x]
    
    visited = {v:0 for v in G.V}
    all_path = list() 
    
    def get_neighbor(G: CD, o:str):
        chs = {(' -> ', ch) for ch in G.ch(o)}
        pas = {(' <- ', pa) for pa in G.pa(o)}
        confs = {(' <-> ', v) for v in G.V if G.is_confounded(v, o)}
        
        return chs | pas | confs
    
    def DFS(G: CD, v: str, y: str, path: list[str]):
        if v == y:
            all_path.append(path + v)
            return
        
        visited[v] = 1
        for d,n in get_neighbor(G, v):
            if not visited[n] :
                DFS(G, n, y, path + v + d)
        visited[v] = 0
    
    DFS(G,x,y,"")
    
    return all_path

def d_connected_path(G:CD, xs: Iterable[str], ys: Iterable[str], zs: Iterable[str] = frozenset()):
    
    xs, ys, zs = wrap(xs), wrap(ys), wrap(zs)
    
    xsys = frozenset(itertools.product(xs,ys))
    
    dps = list()
    for x,y in xsys:
        dps += d_connected_path_helper(G,x,y,zs)
    
    return dps

def d_connected_path_helper(G: CD, x: str, y: str, Zs: Iterable[str] = frozenset()):
    
    Zs = wrap(Zs)
    
    assert x in G.V and y in G.V # x and y are in V
    assert all(z in G.V for z in Zs) # Zs are in V
    assert x not in Zs and y not in Zs
    
    if x == y:
        return [x]
    
    visited = {v:0 for v in G.V}
    visited[x] = 1
    
    def DFS(G: CD, dv: str , v: str, y: str, path: str):
        if v == y:
            trajectoy.append(path)
            return
        
        for dn, n in get_neighbor(G, dv, v, Zs):
            # print(dv,v,__get_neighbor(G, dv, v, Zs))
            if not visited[n] :
                visited[n] = 1
                DFS(G, dn, n, y, path + dn + n)
                visited[n] = 0
    
    trajectoy = list() 
    DFS(G,' <- ',x,y,x)
    
    return trajectoy    

def closure(G: CD, A: Iterable[str], B: Iterable[str], C: Iterable[str]) -> frozenset[str]:
    A,B,C = wrap(A), wrap(B), wrap(C)
    assert not (A & B & C) 
    
    reached = set()
    V = G.An(A|B)
    
    for s in A: 
        if s not in reached: 
            visited = closure_helper(G[V], s, V&C)
            reached |= {v for v in visited if visited[v]}
    return wrap(reached)

def closure_helper(G: CD,
                   x: str,
                   zs: Iterable[str] = frozenset()) -> dict[str,int]:
    zs = wrap(zs)
    
    assert x in G.V and x not in zs 
    assert all(z in G.V for z in zs)
    
    temp = {v:0 for v in G.V}
    visited = {v:0 for v in G.V}
    
    temp[x] = 1
    visited[x] = 1
    
    def DFS(G:CD, dv: str , v: str, path: str):
        if not get_neighbor(G, dv, v, zs) or all(temp[n] for _, n in get_neighbor(G, dv, v, zs)):
            # trajectoy.append(path)
            for v in temp:
                if temp[v]:
                    visited[v] = 1 
            return
        
        for dn, n in get_neighbor(G, dv, v, zs):
            # print(dv,v,__get_neighbor(G, dv, v, zs))
            if not temp[n] :
                temp[n] = 1
                DFS(G, dn, n, path + dn + n)
                temp[n] = 0
    
    # trajectoy = list() 
    DFS(G,' <- ',x,x)
    return visited

def get_neighbor(G:CD, d:str, x:str, Zs:Iterable[str]):
        
        blocked = x in Zs
        colliderables = set_union(G.An(z) for z in Zs)
        
        chs = {(' -> ', ch) for ch in G.ch(x)}
        pas = {(' <- ', pa) for pa in G.pa(x)}
        confs = {(' <-> ', v) for v in G.V if G.is_confounded(v, x)}
        
        nexts = set()
        
        if d == ' -> ' or d == ' <-> ':
            if not blocked:  
                nexts |= chs
                if x in colliderables:
                    nexts |= pas
                    nexts |= confs 
            else:  
                nexts |= pas
                nexts |= confs  
        elif d == ' <- ' :
            if not blocked:  
                nexts |= pas
                nexts |= chs
                nexts |= confs
        else :
            ValueError("Invalid direction")
        
        return nexts

def pcp(G: CD, xs: Iterable[str], ys: Iterable[str]) -> list[str]:
    xs, ys = wrap(xs), wrap(ys)
    def pcp_helper(G: CD, c: str, y: str, X: frozenset[str], path: list[str], paths: list[list[str]]) :
        path = path + [c]
        if c==y:
            paths = paths.append(path)
        
        for n in G.ch(c):
            if n not in path and n not in X:
                pcp_helper(G,n,y,X,path,paths)

    pcps = list()
    xsys = frozenset(itertools.product(xs,ys))
    
    for x,y in xsys:
        pcp_helper(G,x,y,xs,list(),paths:=list())
        pcps += paths[:]
    
    assert pcpSet(G,xs,ys) == wrap(set_union(pcps) - xs)
    return pcps

def pcpSet(G: CD, xs: Iterable[str], ys: Iterable[str]) -> frozenset[str]:        
    xs, ys = wrap(xs), wrap(ys)
    W = (G.do(xs).de(xs) - xs) & G.underbar(xs).An(ys)
    return W

def dpcpSet(G: CD, xs: Iterable[str], ys: Iterable[str]):
    xs, ys = wrap(xs), wrap(ys)
    pcps = pcp(G, xs, ys)
    W = wrap(set_union([p for p in pcps if p[0] in xs and p[-1] in ys]) - xs) # {pcps} - xs
    assert W == pcpSet(G,xs,ys)
    return G.De(W)

def pbd(G:CD, xs: Iterable[str], ys: Iterable[str]) -> CD:
    xs, ys = wrap(xs), wrap(ys)
    assert not (xs & ys) and (xs | ys <= G.V)
    
    G_pbd = deepcopy(G)
    pcps = pcp(G,xs,ys)

    for pcp_ in pcps:
        s,n, e = pcp_[0], pcp_[1], pcp_[-1]  # start, next, end
        if s in xs and e in ys and (s,n) in G_pbd.edges:
            G_pbd = G_pbd.edges_removed([[s,n]])
    return G_pbd

def __cp(G:CD, x:str,y:str) -> list[list[str]]:
    assert isinstance(x,str) and isinstance(y,str)
    
    def cp_helper(G:CD,c:str,y:str,path:list[str],paths:list[list[str]]) :
        path = path + [c]
        if c == y:
            paths = paths.append(path)
        
        for n in G.ch(c):
            if n not in path:
                cp_helper(G,n,y,path,paths)

    cp_helper(G,x,y,list(),paths:=list())
    return paths

def cp(G:CD, xs: Iterable[str], ys: Iterable[str]):
    xs, ys = wrap(xs), wrap(ys)
    xsys = frozenset(itertools.product(xs,ys))
    cps = list()
    for x,y in xsys:
        cps += __cp(G,x,y)
    return cps


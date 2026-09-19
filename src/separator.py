"""Constrained inclusion-minimal separators for causal DAGs and ADMGs.

The observed ADMG is expanded with one latent parent per bidirected edge.
Moralization and greedy deletion are deliberately simple; this module does not
implement the linear-time FindMinSep routine or claim its complexity bound.
"""
from itertools import combinations


def _ancestral_moral_graph(graph, seeds):
    parents = {v: set(graph.pa(v)) for v in graph.V}
    for endpoints in graph.confounded_dict.values():
        latent = object()  # Cannot collide with an observed or user-supplied name.
        parents[latent] = set()
        for v in endpoints:
            parents[v].add(latent)
    ancestors = set(seeds)
    pending = list(ancestors)
    while pending:
        for p in parents[pending.pop()] - ancestors:
            ancestors.add(p)
            pending.append(p)
    adjacency = {v: set() for v in ancestors}
    for v in ancestors:
        ps = parents[v] & ancestors
        for p in ps:
            adjacency[v].add(p)
            adjacency[p].add(v)
        for a, b in combinations(ps, 2):
            adjacency[a].add(b)
            adjacency[b].add(a)
    return adjacency


def _separated(adjacency, left, right, removed):
    reached = set(left)
    pending = list(left)
    while pending:
        v = pending.pop()
        if v in right:
            return False
        for w in adjacency[v] - reached - removed:
            reached.add(w)
            pending.append(w)
    return True


def separated(graph, left, right, conditioned):
    """m-separation, with disjoint observed endpoint and conditioning sets."""
    left, right, conditioned = map(frozenset, (left, right, conditioned))
    if (left & right or left & conditioned or right & conditioned
            or not (left | right | conditioned) <= graph.V):
        raise ValueError("Separation arguments must be disjoint subsets of observed vertices")
    adjacency = _ancestral_moral_graph(graph, left | right | conditioned)
    return _separated(adjacency, left, right, conditioned)


def minimal_separator(graph, left, right, mandatory, candidates):
    """Return mandatory plus an inclusion-minimal supplement, or None.

    Candidates must be ancestors of the endpoints or mandatory vertices. This
    makes the ancestral domain fixed under deletion, so separation is monotone
    in removed vertices in this one moral graph. A single greedy pass therefore
    gives inclusion minimality of the supplement, not minimum cardinality.
    """
    left, right, mandatory, candidates = map(
        frozenset, (left, right, mandatory, candidates))
    if (left & right or (left | right) & (mandatory | candidates)
            or mandatory & candidates
            or not (left | right | mandatory | candidates) <= graph.V):
        raise ValueError("Separator arguments must be disjoint observed vertex sets")
    adjacency = _ancestral_moral_graph(graph, left | right | mandatory)
    if not candidates <= adjacency.keys():
        raise ValueError("Candidates must lie in the mandatory-history ancestor domain")
    selected = set(mandatory | candidates)
    if not _separated(adjacency, left, right, selected):
        return None
    for v in sorted(candidates):
        selected.remove(v)
        if not _separated(adjacency, left, right, selected):
            selected.add(v)
    return frozenset(selected)

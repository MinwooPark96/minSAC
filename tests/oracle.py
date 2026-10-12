"""Small-graph oracle: raw edges and active simple paths, no src helpers."""
from itertools import combinations, product


def reach(vertices, edges, starts, reverse=False):
    found = set(starts)
    while True:
        new = {a if reverse else b for a, b in edges
               if (b if reverse else a) in found} - found
        if not new:
            return found
        found |= new


def separated(vertices, directed, bidirected, x, y, condition):
    ancestors = reach(vertices, directed, condition, reverse=True)
    adjacent = {v: [] for v in vertices}
    for a, b in directed:
        adjacent[a].append((b, False, True))
        adjacent[b].append((a, True, False))
    for a, b in bidirected:
        adjacent[a].append((b, True, True))
        adjacent[b].append((a, True, True))

    def active(v, incoming, seen):
        if v in y:
            return True
        for w, head_here, head_there in adjacent[v]:
            if w in seen:
                continue
            if incoming is not None:
                collider = incoming and head_here
                if (collider and v not in ancestors
                        or not collider and v in condition):
                    continue
            if active(w, head_there, seen | {w}):
                return True
        return False

    return not any(active(v, None, {v}) for v in x)


def partition(vertices, directed, treatments, outcomes):
    parts = [set() for _ in range(len(treatments) + 1)]
    for y in outcomes:
        k = max((i + 1 for i, x in enumerate(treatments)
                 if y in reach(vertices, directed, {x})), default=0)
        parts[k].add(y)
    return parts


def stage(vertices, directed, bidirected, treatments, outcomes, i):
    parts = partition(vertices, directed, treatments, outcomes)
    target = set().union(*parts[i + 1:])
    future = set(treatments[i + 1:])
    d = {(a, b) for a, b in directed if b not in future}
    b = {(a, z) for a, z in bidirected if not ({a, z} & future)}
    # Enumerate directed causal paths, independently of production ancestry code.
    on_path, first_edges = set(), set()

    def visit(path):
        if path[-1] in target:
            on_path.update(path[1:])
            first_edges.add(tuple(path[:2]))
        for a, z in d:
            if a == path[-1] and z not in path:
                visit(path + [z])

    visit([treatments[i]])
    forbidden = reach(vertices, d, on_path)
    return parts, target, d - first_edges, b, forbidden


def valid(vertices, directed, bidirected, treatments, outcomes, blocks):
    if len(blocks) != len(treatments):
        return False
    history = partition(vertices, directed, treatments, outcomes)[0]
    seen = set(treatments) | set(outcomes)
    for i, block in enumerate(blocks):
        block = set(block)
        if not block <= set(vertices) or block & seen:
            return False
        if block & reach(vertices, directed, treatments[i + 1:]):
            return False
        parts, target, d, b, forbidden = stage(
            vertices, directed, bidirected, treatments, outcomes, i)
        if block & forbidden or not separated(
                vertices, d, b, {treatments[i]}, target, history | block):
            return False
        seen |= block
        history |= block | {treatments[i]} | parts[i + 1]
    return True


def allocations(vertices, treatments, outcomes):
    covariates = sorted(set(vertices) - set(treatments) - set(outcomes))
    for assignment in product(range(len(treatments) + 1), repeat=len(covariates)):
        blocks = [set() for _ in treatments]
        for v, k in zip(covariates, assignment):
            if k:
                blocks[k - 1].add(v)
        yield blocks


def deletion_minimal(vertices, directed, bidirected, treatments, outcomes, blocks):
    choices = [[set(s) for n in range(len(block) + 1)
                for s in combinations(sorted(block), n)] for block in blocks]
    return all(not valid(vertices, directed, bidirected, treatments, outcomes, smaller)
               for smaller in product(*choices)
               if tuple(map(frozenset, smaller)) != tuple(map(frozenset, blocks)))


def stagewise_minimal(vertices, directed, bidirected, treatments, outcomes, blocks):
    history = partition(vertices, directed, treatments, outcomes)[0]
    for i, block in enumerate(blocks):
        parts, target, d, b, _ = stage(vertices, directed, bidirected, treatments, outcomes, i)
        for size in range(len(block)):
            for smaller in combinations(sorted(block), size):
                if separated(vertices, d, b, {treatments[i]}, target, history | set(smaller)):
                    return False
        history |= set(block) | {treatments[i]} | parts[i + 1]
    return True

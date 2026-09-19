# Sequential covariate adjustment

Python code accompanying *Complete Graphical Criterion for Sequential Covariate Adjustment in Causal Inference* (NeurIPS 2024), with corrections to the sequential criterion and construction.

The implementation correction is tracked in [issue #4](https://github.com/snu-causality-lab/minSAC/issues/4).

At stage `i`, both the forbidden descendants and the proper back-door graph are computed in `G` after removing incoming arrowheads into **later treatments**, including their incident bidirected edges. Covariates may descend from the current treatment, but cannot descend from later treatments in the original graph. The outcome partition is determined by each outcome's last treatment ancestor in the selected order.

The constructor finalizes each block before computing the next stage's candidates. Its local separator retains mandatory history in the ancestor domain. These changes address the original forbidden-set rejection and the loss of covariates caused by pruning against an outdated history.

## Install and test

The core API requires Python 3.11 or later, NumPy, and NetworkX. From this repository:

```sh
python -m pip install numpy 'networkx>=3.3' pytest
python -m pytest -q
```

`requirements.txt` retains the broader environment used by the original notebooks, including optional R/Jupyter dependencies; it is not needed for the core API or tests.

## Construct a sequence

```python
from src.model import CausalDiagram
from src.sac import SAC

g = CausalDiagram(
    ['C', 'X1', 'X2', 'Y'],
    [('X1', 'X2'), ('C', 'X2'), ('C', 'Y'), ('X2', 'Y')],
)
sac = SAC(g, {'X1', 'X2'}, {'Y'}, order=['X1', 'X2'])
sequence = sac.minsac()
assert sequence == (frozenset(), frozenset({'C'}))
assert sac.criteria()
```

Use string vertex labels. For an ADMG, pass bidirected edges as `(left, right, unique_latent_name)` triples to `CausalDiagram`; directed and bidirected edges may coexist between the same vertices (bows). The directed graph must be acyclic. Do not mutate a graph after constructing `SAC`.

`construct()` returns canonical valid blocks. `minsac()` returns stagewise inclusion-minimal blocks: no strict subset works at a stage with its earlier blocks fixed. It also rules out any valid sequence obtained by componentwise deletion. It does not minimize total cardinality or estimator variance.

Both methods search all observed covariates outside treatments and outcomes. They return a tuple of `frozenset` blocks and update `sac.Z` on success. They return `None` if no admissible sequential-adjustment sequence exists for this order, leaving the previous blocks unchanged. Failure does not mean that the causal effect is unidentifiable by other methods. The return value is explicit; legacy calls that only inspect the mutated `sac.Z` still work after success.

## Check supplied covariates

```python
sac.set_covariates([set(), {'C'}])
assert sac.criteria()
```

`set_covariates(blocks)` installs exact stage assignments atomically. `set_covariate(i, block)` replaces one block, using one-based stage indices. Both reject inadmissible inputs with `ValueError`: unknown vertices, treatment/outcome covariates, repeated covariates across blocks, or descendants of later treatments. An admissible sequence may still fail `criteria()` because it contains forbidden vertices or leaves an active path.

The optional fourth constructor argument `zs` is an **initial covariate union**, allocated at the earliest admissible stage. It is not a restriction on the construction search space. Omit it when constructing a sequence from scratch. Current-treatment descendants are allowed in supplied blocks whenever the stage criterion permits them.

`order=` must contain each treatment exactly once and respect directed ancestry, including paths through other vertices. Without it, a deterministic topological order is selected without enumerating all orders. `set_next()` lazily selects another treatment order and rebuilds the outcome partition and earliest-stage allocation of the current covariate union; previous exact stage assignments are discarded. It returns `False` on exhaustion without changing the current state. `OrderedX.get_all_causal_order()` remains an explicit enumeration API and can have factorial output.

## Verification and computational scope

The tests cover the chain forbidden-set regression, forward reassignment, mandatory-history colliders, impossible bows, atomic failure, supplied-block validation, and treatment order changes. They compare both construction methods against exhaustive covariate allocations on all 64 four-node DAGs in one fixed vertex order and seeded two- and three-treatment ADMG draws, including bows. Every returned minimal sequence is checked for stagewise and componentwise inclusion minimality.

The independent oracle in `tests/oracle.py` uses raw edges, directed-path enumeration, and edge-marked active-simple-path enumeration. It does not call production graph transformations, ancestry, separation, or construction helpers. These are finite graphical regression checks, not a proof of the universal formula identity or a running-time benchmark.

Production separation uses a latent DAG expansion and ancestral moralization. The constrained separator deletes optional vertices greedily in one fixed ancestor domain. This implementation **does not claim linear time**: moralization may densify the graph, greedy deletion repeats reachability checks, and existing proper-causal-path helpers enumerate paths. A linear-time constrained-separator implementation would be a separate performance change.

The graph representation builds on the structural causal bandit library by [Sanghack Lee](https://github.com/sanghack81/SCMMAB-NIPS2018).

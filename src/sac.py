"""Sequential adjustment after interventions on future treatments."""
from typing import Iterable

from .model import CausalDiagram as CD
from .utils import fzset_union, wrap
from .order import OrderedX, OrderedY, Orderedsac
from .basis import SequentialAdjustment
from .ac import AdjustmentCriterion
from . import navigator, separator


class SequentialAdjustmentCriterion(SequentialAdjustment):
    """Check supplied blocks or construct a valid sequence for a fixed order.

    ``zs`` is an initial covariate union, allocated at the earliest admissible
    stage. Use ``set_covariates`` to supply exact blocks. Construction searches
    all observed vertices outside treatments/outcomes, regardless of ``zs``.
    ``order`` explicitly selects a topological treatment order; otherwise one
    deterministic order is used. Graphs are expected to remain unchanged.
    """

    def __init__(self, G: CD, xs: Iterable[str], ys: Iterable[str],
                 zs: Iterable[str] = frozenset(), *, order=None) -> None:
        self.G, self.V = G, G.V
        self.xs, self.ys, self.zs = wrap(xs), wrap(ys), wrap(zs)
        if not self.xs or not (self.xs | self.ys | self.zs) <= G.V:
            raise ValueError("Use known vertices and at least one treatment")
        if self.xs & self.ys or self.xs & self.zs or self.ys & self.zs:
            raise ValueError("Treatments, outcomes, and covariates must be pairwise disjoint")
        self.X = OrderedX(G, self.xs)
        if order is not None:
            self.X.set_partition(order)
        self.m = self.X.m
        self.Y = OrderedY(G, self.X, self.ys)
        self.Z = Orderedsac(G, self.X, self.zs)
        initial_order = tuple(next(iter(self.X[i])) for i in range(1, self.m + 1))
        self.order_iter = (o for o in self.X.iter_causal_orders() if tuple(o) != initial_order)

    def set_next(self) -> bool:
        """Select another order lazily and reallocate the current covariate union.

        Exact previous stage assignments are discarded because stages change.
        Exhaustion returns False and leaves the current sequence unchanged.
        """
        next_order = next(self.order_iter, None)
        if next_order is None:
            return False
        self.X.set_partition(next_order)
        self.Y = OrderedY(self.G, self.X, self.ys)
        self.Z = Orderedsac(self.G, self.X, self.zs)
        return True

    def _stage_index(self, i):
        if not isinstance(i, int) or not 1 <= i <= self.m:
            raise ValueError("Stage must lie between 1 and the number of treatments")

    def _validate_blocks(self, blocks):
        blocks = tuple(wrap(block) for block in blocks)
        if len(blocks) != self.m:
            raise ValueError("Supply exactly one covariate block per treatment")
        seen = self.xs | self.ys
        for i, block in enumerate(blocks, 1):
            if not block <= self.V or block & seen:
                raise ValueError("Covariate blocks must be disjoint observed nontreatment, nonoutcome sets")
            if block & self.G.De(self.X.geq(i + 1)):
                raise ValueError("Covariates cannot descend from a later treatment")
            seen |= block
        return blocks

    def set_covariates(self, blocks):
        """Atomically install an admissible sequence; validity is checked separately."""
        blocks = self._validate_blocks(blocks)
        union = fzset_union(blocks)
        new_z = Orderedsac(self.G, self.X, union)
        new_z.P = list(blocks)
        self.Z, self.zs = new_z, union

    def set_covariate(self, i: int, zi: Iterable[str] = frozenset()):
        self._stage_index(i)
        blocks = list(self.Z.P)
        blocks[i - 1] = wrap(zi)
        self.set_covariates(blocks)

    def dpcpSet(self, i: int):
        """Forbidden descendants, entirely in the future-intervened graph."""
        self._stage_index(i)
        future = self.G.do(self.X.geq(i + 1))
        return navigator.dpcpSet(future, self.X[i], self.Y.geq(i))

    def condition1(self, i: int) -> bool:
        return not self.Z[i] & self.dpcpSet(i)

    def condition2(self, i: int) -> bool:
        self._stage_index(i)
        return separator.separated(self.psbd(i), self.X[i], self.Y.geq(i),
                                   self.H(i - 1) | self.Z[i])

    def criteria(self) -> bool:
        try:
            self._validate_blocks(self.Z.P)
        except ValueError:
            return False
        return all(self.condition1(i) and self.condition2(i)
                   for i in range(1, self.m + 1))

    def psbd(self, i: int):
        self._stage_index(i)
        return navigator.pbd(self.G.do(self.X.geq(i + 1)), self.X[i], self.Y.geq(i))

    def F(self, i: int):
        self._stage_index(i)
        return (self.xs | self.ys | self.H(i - 1) | self.dpcpSet(i)
                | self.G.De(self.X.geq(i + 1)))

    def _construct(self, minimal):
        history = self.Y[0]
        blocks = []
        for i in range(1, self.m + 1):
            target = self.Y.geq(i)
            if not target:
                block = frozenset()
            else:
                stage_graph = self.psbd(i)
                forbidden = (self.xs | self.ys | history | self.dpcpSet(i)
                             | self.G.De(self.X.geq(i + 1)))
                candidates = stage_graph.An(self.X[i] | target | history) - forbidden
                if minimal:
                    selected = separator.minimal_separator(
                        stage_graph, self.X[i], target, history, candidates)
                    if selected is None:
                        return None
                    block = selected - history
                else:
                    if not separator.separated(stage_graph, self.X[i], target,
                                               history | candidates):
                        return None
                    block = candidates
            blocks.append(block)
            history |= self.X[i] | self.Y[i] | block
        self.set_covariates(blocks)
        return tuple(blocks)

    def construct(self):
        """Install/return canonical valid blocks, or None without changing state."""
        return self._construct(minimal=False)

    def minsac(self):
        """Install/return inclusion-minimal blocks, or None without changing state.

        Earlier blocks are finalized before later candidates are computed.
        Minimality is by set inclusion, not cardinality or estimator variance.
        The implementation does not claim linear time per treatment.
        """
        return self._construct(minimal=True)

    def to_ac(self) -> AdjustmentCriterion:
        return AdjustmentCriterion(G=self.G, xs=self.xs, ys=self.ys,
                                   zs=fzset_union(self.Z.P))

    def __getitem__(self, i: int):
        return self.psbd(i)

    def __repr__(self):
        return f"SequentialAdjustmentCriterion({self.G},\n{self.X},\n{self.Y},\n{self.Z})"

    __str__ = __repr__


SAC = SequentialAdjustmentCriterion

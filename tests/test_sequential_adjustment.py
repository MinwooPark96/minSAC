from itertools import combinations
import random

import pytest

from src.model import CausalDiagram
from src.order import OrderedX
from src.sac import SAC
from tests.oracle import allocations, deletion_minimal, separated, stagewise_minimal, valid


def graph(vertices, directed=(), bidirected=()):
    return CausalDiagram(vertices, directed,
                         [(a, b, f"__u{i}") for i, (a, b) in enumerate(bidirected)])


def blocks(sac):
    return tuple(sac.Z[i] for i in range(1, sac.m + 1))


def test_chain_forbidden_uses_future_intervention():
    g = graph(["X1", "Z", "X2", "Y"],
              [("X1", "Z"), ("Z", "X2"), ("X2", "Y")])
    sac = SAC(g, ["X1", "X2"], {"Y"}, {"Z"})
    sac.set_covariate(2, set())
    sac.set_covariate(1, {"Z"})
    assert "Z" not in sac.dpcpSet(1)
    assert sac.criteria()


def test_pruned_covariate_is_reconsidered_later():
    edges = [("X1", "X2"), ("C", "X2"), ("C", "Y"), ("X2", "Y")]
    sac = SAC(graph(["C", "X1", "X2", "Y"], edges), ["X1", "X2"], {"Y"}, {"C"})
    result = sac.minsac()
    assert blocks(sac) == (frozenset(), frozenset({"C"}))
    assert result == blocks(sac)
    assert sac.criteria()


def test_mandatory_history_collider_is_in_separator_domain():
    vertices = ["A", "B", "C", "X", "Y"]
    edges = [("A", "X"), ("A", "C"), ("B", "C"), ("B", "Y"), ("X", "Y")]
    sac = SAC(graph(vertices, edges), {"X"}, {"C", "Y"}, {"A", "B"})
    sac.minsac()
    assert valid(vertices, edges, [], ["X"], {"C", "Y"}, blocks(sac))
    assert sum(map(len, blocks(sac))) == 1


def test_construction_uses_all_observed_covariates_by_default():
    sac = SAC(graph(["C", "X", "Y"], [("C", "X"), ("C", "Y"), ("X", "Y")]),
              {"X"}, {"Y"})
    assert sac.minsac() == (frozenset({"C"}),)
    assert sac.criteria()


@pytest.mark.parametrize("method", ["construct", "minsac"])
def test_impossible_bow_reports_failure_without_mutation(method):
    sac = SAC(graph(["Z", "X", "Y"], [("X", "Y")], [("X", "Y")]),
              {"X"}, {"Y"}, {"Z"})
    original = blocks(sac)
    assert getattr(sac, method)() is None
    assert blocks(sac) == original
    assert not sac.criteria()


def test_explicit_order_and_supplied_blocks():
    g = graph(["X1", "X2", "Y", "C"], [("X1", "Y"), ("X2", "C")])
    sac = SAC(g, {"X1", "X2"}, {"Y"}, order=["X2", "X1"])
    assert [next(iter(sac.X[i])) for i in [1, 2]] == ["X2", "X1"]
    sac.set_covariates([{"C"}, set()])
    assert blocks(sac) == (frozenset({"C"}), frozenset())
    assert sac.criteria()
    original = blocks(sac)
    with pytest.raises(ValueError):
        sac.set_covariates([{"C"}, {"C"}])
    assert blocks(sac) == original


def test_topological_order_validation_uses_paths_through_other_vertices():
    g = graph(["X1", "M", "X2", "Y"], [("X1", "M"), ("M", "X2")])
    for order in [["X2", "X1"], ["X1", "X1"], ["X1"]]:
        with pytest.raises(ValueError):
            SAC(g, {"X1", "X2"}, {"Y"}, order=order)
    ordered = OrderedX(g, {"X1", "X2"})
    with pytest.raises(ValueError):
        ordered.set_partition(["X2", "X1"])
    assert ordered & {"X1"} == frozenset({"X1"})


@pytest.mark.parametrize("xs,ys,zs", [({"X"}, {"X"}, set()),
                                     ({"X"}, {"Y"}, {"X"}),
                                     ({"X"}, {"Y"}, {"unknown"}),
                                     (set(), {"Y"}, set())])
def test_invalid_vertex_inputs_raise(xs, ys, zs):
    with pytest.raises(ValueError):
        SAC(graph(["X", "Y"]), xs, ys, zs)


def test_supplied_future_descendants_and_unknown_vertices_rejected():
    sac = SAC(graph(["X1", "X2", "C", "Y"], [("X2", "C")]),
              {"X1", "X2"}, {"Y"}, order=["X1", "X2"])
    for bad in [[{"C"}, set()], [{"unknown"}, set()], [{"Y"}, set()], [set()]]:
        with pytest.raises(ValueError):
            sac.set_covariates(bad)


def test_next_order_rebuilds_outcome_and_covariate_partitions():
    g = graph(["X1", "X2", "Y", "C"], [("X1", "Y"), ("X2", "C")])
    sac = SAC(g, {"X1", "X2"}, {"Y"}, {"C"})
    old_order = tuple(next(iter(sac.X[i])) for i in [1, 2])
    assert sac.set_next()
    new_order = tuple(next(iter(sac.X[i])) for i in [1, 2])
    assert new_order != old_order
    assert sac.Y[new_order.index("X1") + 1] == frozenset({"Y"})
    assert sac.Z[new_order.index("X2") + 1] == frozenset({"C"})
    assert not sac.set_next()


def test_empty_outcome_stage_is_vacuous_and_returns_empty_block():
    sac = SAC(graph(["X1", "X2", "Y"], [("X1", "Y")]),
              {"X1", "X2"}, {"Y"}, order=["X1", "X2"])
    assert sac.minsac() == (frozenset(), frozenset())
    assert sac.criteria()


def test_constructor_does_not_enumerate_all_orders(monkeypatch):
    def forbidden_enumeration(*args, **kwargs):
        raise AssertionError("Initialization must not enumerate all orders")
    monkeypatch.setattr(OrderedX, "get_all_causal_order", forbidden_enumeration)
    sac = SAC(graph(["X1", "X2", "Y"]), {"X1", "X2"}, {"Y"})
    assert sac.minsac() == (frozenset(), frozenset())


def test_cycles_and_invalid_stage_indices_raise():
    with pytest.raises(ValueError):
        SAC(graph(["X", "Y"], [("X", "Y"), ("Y", "X")]), {"X"}, {"Y"})
    sac = SAC(graph(["X", "Y"]), {"X"}, {"Y"})
    for i in [0, 2, -1]:
        with pytest.raises(ValueError):
            sac.set_covariate(i, set())


@pytest.mark.parametrize("method", ["construct", "minsac"])
def test_failure_at_later_stage_is_atomic(method):
    sac = SAC(graph(["Z", "X1", "X2", "Y"], [("X1", "X2"), ("X2", "Y")],
                    [("X2", "Y")]), {"X1", "X2"}, {"Y"}, {"Z"})
    old_blocks, old_union = blocks(sac), sac.zs
    assert getattr(sac, method)() is None
    assert blocks(sac) == old_blocks and sac.zs == old_union


def test_mandatory_history_separator_and_bow_separation():
    from src.separator import minimal_separator, separated as production_separated
    vertices = ["A", "B", "C", "X", "Y"]
    edges = [("A", "X"), ("A", "C"), ("B", "C"), ("B", "Y")]
    g = graph(vertices, edges)
    result = minimal_separator(g, {"X"}, {"Y"}, {"C"}, {"A", "B"})
    assert result is not None and "C" in result and len(result) == 2
    assert separated(vertices, edges, [], {"X"}, {"Y"}, result)
    with pytest.raises(ValueError):
        minimal_separator(graph(["X", "Y", "Z"]), {"X"}, {"Y"}, set(), {"Z"})
    for conditioning in [set(), {"C"}]:
        bidirected = [("X", "Y"), ("A", "X")]
        assert production_separated(graph(vertices, edges, bidirected), {"X"}, {"Y"}, conditioning) == separated(
            vertices, edges, bidirected, {"X"}, {"Y"}, conditioning)


def check_collection(vertices, treatments, outcomes, cases):
    for directed, bidirected in cases:
        g = graph(vertices, directed, bidirected)
        sac = SAC(g, treatments, outcomes, order=treatments)
        solutions = []
        for allocation in allocations(vertices, treatments, outcomes):
            expected = valid(vertices, directed, bidirected, treatments, outcomes, allocation)
            try:
                sac.set_covariates(allocation)
            except ValueError:
                assert not expected
            else:
                assert sac.criteria() == expected, (directed, bidirected, allocation)
            if expected:
                solutions.append(allocation)
        for method in ["construct", "minsac"]:
            result = getattr(sac, method)()
            assert (result is not None) == bool(solutions), (directed, bidirected, method)
            if result is not None:
                assert valid(vertices, directed, bidirected, treatments, outcomes, result)
                if method == "minsac":
                    assert deletion_minimal(vertices, directed, bidirected, treatments, outcomes, result)
                    assert stagewise_minimal(vertices, directed, bidirected, treatments, outcomes, result)


def test_exhaustive_four_node_dags_against_independent_oracle():
    vertices = ["C", "X1", "X2", "Y"]
    possible = list(combinations(vertices, 2))
    cases = [([edge for j, edge in enumerate(possible) if mask & (1 << j)], [])
             for mask in range(1 << len(possible))]
    check_collection(vertices, ["X1", "X2"], {"Y"}, cases)


def test_seeded_admgs_including_bows_against_independent_oracle():
    vertices = ["A", "X1", "B", "C", "X2", "Y"]
    possible = list(combinations(vertices, 2))
    rng = random.Random(20260919)
    cases = [([edge for edge in possible if rng.random() < .35],
              [edge for edge in possible if rng.random() < .18]) for _ in range(128)]
    assert any(set(d) & set(b) for d, b in cases)
    check_collection(vertices, ["X1", "X2"], {"C", "Y"}, cases)


def test_three_treatment_admgs_against_independent_oracle():
    vertices = ["A", "X1", "X2", "C", "X3", "Y"]
    possible = list(combinations(vertices, 2))
    rng = random.Random(20260920)
    cases = [([edge for edge in possible if rng.random() < .35],
              [edge for edge in possible if rng.random() < .18]) for _ in range(64)]
    check_collection(vertices, ["X1", "X2", "X3"], {"C", "Y"}, cases)

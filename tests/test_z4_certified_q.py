"""A finite-q compression witness must not silently become an all-q scheme."""
import pytest


@pytest.mark.parametrize('certified,requested', [(8, 7), (True, 1), (8.0, 8), (0, 0), ('8', 8), (8, 8.0), ([], 8), ([2, 8], 7)])
def test_wrong_or_malformed_q_is_rejected_before_building_tables(monkeypatch, certified, requested):
    from cwschemes._engine import z4_quotient as z4
    monkeypatch.setattr(z4, 'relations_object', lambda path: ({'certified_q': certified}, [], 'unused'))
    with pytest.raises(ValueError, match='certified only'):
        z4.replay_document('finite-q', requested)
    with pytest.raises(ValueError, match='certified only'):
        z4.JointPacket(None, requested, 'finite-q')


def test_matching_q_and_legacy_replays_keep_their_scope(monkeypatch):
    from cwschemes._engine import z4_quotient as z4
    d = dict(m=3, rows=[], initial_pair=[0, 1], steps=[dict(
        pair=[0, 2], kept=[0], deleted=[], zero_patterns=[], blocks={})])
    monkeypatch.setattr(z4, 'relations_object', lambda path: (d, [], 'unused'))
    for q in (2, 8, 17):
        assert z4.replay_document('legacy', q)['q'] == q
    d['certified_q'] = 8
    assert z4.replay_document('finite-q', 8)['complete']
    d['certified_q'] = list(range(2, 26))
    for q in (2, 8, 17, 25):
        assert z4.replay_document('historical-replays', q)['q'] == q


def test_frozen_pattern_counts_are_separate_from_coefficient_scope():
    from cwschemes._engine.z4_checks import closed_form_rank
    polys = dict(even=[1,17,122,392,493,127,2,2,0,0],
                 odd=[1,17,122,392,493,127,-2,-2,0,0])
    assert closed_form_rank(8, 'odd', polys) == 198679264
    assert closed_form_rank(2, 'odd', polys) == 15838
    assert closed_form_rank(8, 'even', polys) == 198679840
    assert closed_form_rank(8, 'odd') == 198683360
    with pytest.raises(AssertionError):
        closed_form_rank(2, 'odd', {'odd': polys['odd'][:-1] + [1]})

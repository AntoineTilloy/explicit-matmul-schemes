"""The compact current-LITA accessor agrees with the pinned serializer."""
from fractions import Fraction
import hashlib
from pathlib import Path
import pytest

from cwschemes._engine.current_lita import CurrentLITA
from cwschemes._engine.sources import lita_current as L


@pytest.mark.parametrize('n',[8,14])
def test_all_current_lita_terms(n):
    digest=hashlib.sha256(Path(L.__file__).read_bytes()).hexdigest()
    compact=CurrentLITA(n,digest);expected=L.lita(n)
    assert compact.rank==expected.rank==L.lita_rank(n)
    arrays={}
    for name in 'UVW':arrays.update(getattr(expected,name).arrays(name.lower()))
    for t in range(compact.rank):
        forms=compact.term(t)
        for mode,axis in enumerate('uvw'):
            ptr=arrays[axis+'_indptr'];want={}
            for j in range(int(ptr[t]),int(ptr[t+1])):
                coord=divmod(int(arrays[axis+'_indices'][j]),n)
                want[coord]=Fraction(int(arrays[axis+'_numerators'][j]),int(arrays[axis+'_denominators'][j]))
            assert forms[mode]==want,(n,t,axis)


def test_current_lita_source_binding():
    with pytest.raises((AssertionError,ValueError)):
        CurrentLITA(8,'0'*64)

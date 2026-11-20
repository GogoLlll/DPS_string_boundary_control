import pytest

from wavesys.params import StringParams

def test_travel_time():
    assert StringParams(1.5, 2.0).t0 == pytest.approx(0.75)

@pytest.mark.parametrize("l, a", [(0, 1), (-1, 1), (1, 0), (1, -2), (float("nan"), 1)])
def test_rejects_non_positive(l, a):
    with pytest.raises(ValueError):
        StringParams(l, a)

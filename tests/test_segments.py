import numpy as np

from wavesys.segments import before, split_time

def test_split_reconstructs_argument():
    t0 = 0.75
    z = np.linspace(0, 10, 1001)
    k, s = split_time(z, t0)
    assert np.all((s >= 0) & (s < t0))
    assert np.allclose(k * t0 + s, z, atol=1e-11)

def test_joint_takes_right_side():
    t0 = 0.08                                     # 0.4 / 5: кратные t0 не представимы точно
    for m in range(1, 60):
        k, s = split_time(m * t0, t0)
        assert (int(k), float(s)) == (m, 0.0)
        k, s = split_time(m * t0 * (1 - 1e-15), t0)
        assert (int(k), float(s)) == (m, 0.0)

def test_point_clearly_left_of_joint_stays_left():
    t0 = 0.5
    k, s = split_time(1.0 - 1e-6, t0)
    assert int(k) == 1 and abs(float(s) - (0.5 - 1e-6)) < 1e-12

def test_keeps_shape():
    k, s = split_time(np.zeros((3, 4)) + 0.2, 0.5)
    assert k.shape == s.shape == (3, 4)

def test_before_is_exact_ordering():
    assert before(1, 0.3, 2, 0.0)
    assert before(2, 0.1, 2, 0.2)
    assert not before(2, 0.2, 2, 0.2)
    assert not before(3, 0.0, 2, 0.9)

import numpy as np
import pytest

from wavesys.functions import Expr, ExpressionError, Primitive

def test_arithmetic_and_power():
    f = Expr("2*x^2 - x/4 + 1", "x")
    x = np.linspace(-1, 2, 7)
    assert np.allclose(f(x), 2 * x ** 2 - x / 4 + 1)

def test_constants_and_parameters():
    f = Expr("sin(pi*x/l) + t0", "x", {"l": 2.0, "t0": 0.5})
    assert f(1.0) == pytest.approx(1.5)

def test_constant_formula_returns_array_of_input_shape():
    f = Expr("0.5", "t")
    assert f(np.zeros((2, 3))).shape == (2, 3)
    assert np.all(f(np.zeros(4)) == 0.5)

def test_helper_functions():
    x = np.array([-1.0, 0.0, 0.5, 1.0, 2.0])
    assert np.array_equal(Expr("pos(x)", "x")(x), [0, 0, 0.5, 1, 2])
    assert np.array_equal(Expr("ind(x, 0, 1)", "x")(x), [0, 1, 1, 0, 0])
    assert np.array_equal(Expr("step(x)", "x")(x), [0, 1, 1, 1, 1])

@pytest.mark.parametrize("text", [
    "", "   ", "x +", "y + 1", "foo(x)", "__import__('os')", "x.real", "x[0]", "'abc'", "lambda: x",
    "sin(x, key=1)", "(x for x in x)", "True",
])
def test_rejects_bad_formulas(text):
    with pytest.raises(ExpressionError):
        Expr(text, "x")

def test_error_message_names_the_unknown_symbol():
    with pytest.raises(ExpressionError, match="zeta"):
        Expr("x + zeta", "x")

def test_primitive_of_smooth_function():
    F = Primitive(lambda x: np.cos(3 * x), 0.0, 2.0)
    x = np.linspace(0, 2, 97)
    assert np.allclose(F(x), np.sin(3 * x) / 3, atol=1e-13)
    assert F(0.0) == 0.0

def test_primitive_of_function_with_jump():
    F = Primitive(lambda x: np.where(x < 0.3, 1.0, -2.0), 0.0, 1.0)
    x = np.linspace(0, 1, 51)
    exact = np.where(x < 0.3, x, 0.3 - 2 * (x - 0.3))
    assert np.allclose(F(x), exact, atol=1e-10)         # ячейка со скачком измельчается автоматически

def test_primitive_keeps_shape_and_clips():
    F = Primitive(lambda x: np.ones(np.shape(x)), 0.0, 1.0)
    assert F(np.zeros((2, 5))).shape == (2, 5)
    assert F(5.0) == pytest.approx(1.0)

def test_primitive_of_function_with_kink_from_formula():
    f = Expr("-0.6*pi/1.2*sin(pi*x/0.6)*ind(x, 0, 0.6)", "x")
    F = Primitive(f, 0.0, 1.0)
    x = np.linspace(0, 1, 2001)
    exact = np.where(x < 0.6, 0.6 * (np.cos(np.pi * x / 1.2) ** 2 - 1), -0.6)
    assert np.allclose(F(x), exact, atol=1e-12)

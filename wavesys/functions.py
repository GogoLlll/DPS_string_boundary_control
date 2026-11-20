import ast

import numpy as np

class ExpressionError(ValueError):
    """Ошибка в формуле"""

def _pos(x):
    return np.maximum(x, 0.0)

def _ind(x, lo, hi):
    return np.where((x >= lo) & (x < hi), 1.0, 0.0)

def _step(x):
    return np.where(x >= 0, 1.0, 0.0)

ALLOWED_FUNCTIONS = {
    "sin": np.sin, "cos": np.cos, "tan": np.tan, "exp": np.exp, "log": np.log, "sqrt": np.sqrt,
    "abs": np.abs, "sinh": np.sinh, "cosh": np.cosh, "tanh": np.tanh, "sign": np.sign,
    "min": np.minimum, "max": np.maximum, "where": np.where,
    "pos": _pos,      # pos(x) = max(x, 0)
    "ind": _ind,      # ind(x, lo, hi) = 1 при lo <= x < hi, иначе 0
    "step": _step,    # step(x) = 1 при x >= 0, иначе 0
}
ALLOWED_CONSTANTS = {"pi": np.pi, "e": np.e}

_ALLOWED_NODES = (
    ast.Expression, ast.BinOp, ast.UnaryOp, ast.Call, ast.Name, ast.Constant, ast.Compare, ast.Load,
    ast.Add, ast.Sub, ast.Mult, ast.Div, ast.Pow, ast.Mod, ast.USub, ast.UAdd,
    ast.Lt, ast.Gt, ast.LtE, ast.GtE,
)

class Expr:
    def __init__(self, text, var, constants=None):
        self.text = str(text).strip()
        self.var = var
        self.constants = dict(ALLOWED_CONSTANTS)
        self.constants.update(constants or {})
        if not self.text:
            raise ExpressionError("Пустая формула")
        source = self.text.replace("^", "**")
        try:
            tree = ast.parse(source, mode="eval")
        except SyntaxError as exc:
            raise ExpressionError(f"Не удалось разобрать формулу «{self.text}»: {exc.msg}") from None
        for node in ast.walk(tree):
            if not isinstance(node, _ALLOWED_NODES):
                raise ExpressionError(f"В формуле «{self.text}» недопустимая конструкция: {type(node).__name__}")
            if isinstance(node, ast.Call):
                if not isinstance(node.func, ast.Name) or node.func.id not in ALLOWED_FUNCTIONS or node.keywords:
                    name = getattr(node.func, "id", "?")
                    raise ExpressionError(f"В формуле «{self.text}» неизвестная функция: {name}")
            elif isinstance(node, ast.Name):
                if node.id != var and node.id not in ALLOWED_FUNCTIONS and node.id not in self.constants:
                    raise ExpressionError(
                        f"В формуле «{self.text}» неизвестное имя: {node.id}. "
                        f"Переменная здесь — {var}; параметры: {', '.join(sorted(self.constants))}")
            elif isinstance(node, ast.Constant):
                if not isinstance(node.value, (int, float)) or isinstance(node.value, bool):
                    raise ExpressionError(f"В формуле «{self.text}» допустимы только числа")
        self._code = compile(tree, "<формула>", "eval")
        try:
            self(np.array([0.5]))
        except ExpressionError:
            raise
        except Exception as exc:
            raise ExpressionError(f"Формула «{self.text}» не вычисляется: {exc}") from None

    def __call__(self, x):
        x = np.asarray(x, dtype=float)
        scope = dict(ALLOWED_FUNCTIONS)
        scope.update(self.constants)
        scope[self.var] = x
        with np.errstate(all="ignore"):
            value = eval(self._code, {"__builtins__": {}}, scope)
        return np.zeros(x.shape) + np.asarray(value, dtype=float)

    def __repr__(self):
        return f"Expr({self.text!r}, {self.var!r})"

class Primitive:
    def __init__(self, f, lo, hi, cells=1024, depth=44, tol=1e-13):
        self.f, self.lo, self.hi = f, float(lo), float(hi)
        self._nodes, self._weights = np.polynomial.legendre.leggauss(8)
        edges = np.linspace(self.lo, self.hi, int(cells) + 1)
        left, right = edges[:-1], edges[1:]
        whole = self._cell_integrals(left, right)
        scale = max(float(np.sum(np.abs(whole))), 1e-300)
        done_l, done_r, done_i = [], [], []
        for _ in range(depth):
            mid = 0.5 * (left + right)
            a, b = self._cell_integrals(left, mid), self._cell_integrals(mid, right)
            rough = np.abs(whole - (a + b)) > tol * scale
            done_l.append(left[~rough]); done_r.append(right[~rough]); done_i.append((a + b)[~rough])
            if not np.any(rough):
                left = right = whole = np.array([])
                break
            left, right = np.concatenate([left[rough], mid[rough]]), np.concatenate([mid[rough], right[rough]])
            whole = np.concatenate([a[rough], b[rough]])
        done_l.append(left); done_r.append(right); done_i.append(whole)
        left, integrals = np.concatenate(done_l), np.concatenate(done_i)
        order = np.argsort(left)
        self._edges = left[order]
        self._cum = np.concatenate([[0.0], np.cumsum(integrals[order])])
        self.cells = len(self._edges)

    def _cell_integrals(self, left, right):
        half = (right - left) / 2
        pts = (left + half)[:, None] + half[:, None] * self._nodes[None, :]
        return (np.asarray(self.f(pts), dtype=float) * self._weights[None, :]).sum(axis=1) * half

    def __call__(self, x):
        x = np.asarray(x, dtype=float)
        xc = np.clip(x, self.lo, self.hi)
        idx = np.clip(np.searchsorted(self._edges, xc, side="right") - 1, 0, self.cells - 1)
        left = self._edges[idx]
        half = (xc - left) / 2
        pts = (left + half)[..., None] + half[..., None] * self._nodes
        local = (np.asarray(self.f(pts), dtype=float) * self._weights).sum(axis=-1) * half
        return self._cum[idx] + local

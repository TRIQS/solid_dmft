# https://stackoverflow.com/a/30516254

import ast
import math


class MathExpr(object):
    allowed_nodes = (
        ast.Module,
        ast.Expr,
        ast.Load,
        ast.Expression,
        ast.Add,
        ast.Sub,
        ast.UnaryOp,
        ast.Num,
        ast.BinOp,
        ast.Mult,
        ast.Div,
        ast.Pow,
        ast.BitOr,
        ast.BitAnd,
        ast.BitXor,
        ast.USub,
        ast.UAdd,
        ast.FloorDiv,
        ast.Mod,
        ast.LShift,
        ast.RShift,
        ast.Invert,
        ast.Call,
        ast.Name,
    )

    functions = {
        "abs": abs,
        "complex": complex,
        "min": min,
        "max": max,
        "pow": pow,
        "round": round,
    }
    functions.update(
        {key: value for (key, value) in vars(math).items() if not key.startswith("_")}
    )

    def __init__(self, expr):
        if any(elem in expr for elem in "\n#"):
            raise ValueError(expr)

        node = ast.parse(expr.strip(), mode="eval")
        for curr in ast.walk(node):
            if not isinstance(curr, self.allowed_nodes):
                raise ValueError(curr)

        self.code = compile(node, "<string>", "eval")

    def __call__(self, **kwargs):
        return eval(self.code, {"__builtins__": None}, {**self.functions, **kwargs})

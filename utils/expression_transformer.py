import sympy
from sympy import Piecewise, piecewise_fold, integrate, S, zoo


#Evaluates definite integral F(x) - F(a) rigorously to prevent Constant Drift.

def get_piecewise_integral(expr, x, a):
    
    if expr == 0:
        return S.Zero
    
    F = integrate(expr, x).rewrite(Piecewise)
    try:
        val_at_a = F.subs(x, a)
        if not val_at_a.has(zoo, S.NaN, S.Infinity, sympy.oo, -sympy.oo):
            from sympy import nsimplify
            F = nsimplify((F - val_at_a).rewrite(Piecewise)).expand()
    except (ValueError, TypeError):
        pass
    return F

# Recursively counts total branches in a Piecewise expression.

def _count_piecewise_args(expr):
    
    if not isinstance(expr, Piecewise):
        return 1
    count = 0
    for e, _ in expr.args:
        count += _count_piecewise_args(e)
    return count

#Extracts the Piecewise rewrite and integration logic for transcendental functions.

def convert_to_piecewise_integral(f, x, boundary):
    
    try:
        pw = piecewise_fold(f.rewrite(Piecewise))

        if _count_piecewise_args(pw) > 50:
            from sympy import Integral
            return Integral(f, x)

        if isinstance(pw, Piecewise):
            new_args = []
            for expr, cond in pw.args:
                F_adjusted = get_piecewise_integral(expr, x, boundary)
                new_args.append((F_adjusted, cond))

            res = Piecewise(*new_args)
            return piecewise_fold(res)
    except (RecursionError, ValueError):
        pass
    
    return None

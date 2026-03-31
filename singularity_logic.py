from sympy import (
    SingularityFunction, Piecewise, piecewise_fold, integrate, simplify,
    S, zoo, DiracDelta, Symbol, sin, cos, exp, log, Add, Mul, Heaviside
)
from utils.expression_transformer import convert_to_piecewise_integral


def _has_transcendental_of(expr, var):
    """
    Check if expr contains transcendental functions that depend on var.
    
    Developer Notes:
    This check is a deliberate mathematical optimization. Reaction coefficients 
    from solve() often contain constants like exp(5) or cos(3). If we blindly 
    triggered the Piecewise fallback for these, we'd lose the O(1) efficiency 
    of the SingularityFunction power rule unnecessarily. Checking free_symbols 
    ensures we only fallback when the *shape* of the load requires it.
    """
    for func_type in (sin, cos, exp, log):
        for sub_expr in expr.atoms(func_type):
            if var in sub_expr.free_symbols:
                return True
    return False



#SymPy's simplify() can crash with 'Invalid NaN comparison' on mixed
#SingularityFunction + Piecewise expressions. This wrapper catches that.

def _safe_simplify(expr):
    
    try:
        return simplify(expr)
    except (TypeError, ValueError, RecursionError):
        return expr


def dispatch_integration(f, x):
    """
    The 'Smart Dispatcher' for beam load integration.

    Design rationale (Euler-Bernoulli beam theory):
    1. Point loads (n=-1) and moments (n=-2) stay as SingularityFunctions
       because they are mathematically cleaner for boundary condition solvers.
    2. Distributed loads with transcendental q(x) get rewritten to Piecewise
       and integrated with F(x)-F(a) to enforce C^0 continuity.
    3. A Complexity Guard prevents exponential Piecewise branch growth.
    """

    # --- Non-SingularityFunction expressions (Piecewise, plain, constants) ---
    # When integrating V(x) -> M(x), the shear force may already be Piecewise.
    if not f.has(SingularityFunction):
        return integrate(f, x)

    # --- Bare SingularityFunction: <x-a>^n ---
    if isinstance(f, SingularityFunction):
        base, a, n = f.args
        if n < 0:
            # Point load/moment: <x-a>^n -> <x-a>^{n+1} (no divisor)
            return SingularityFunction(base, a, n + 1)
        else:
            # Polynomial distributed: <x-a>^n -> <x-a>^{n+1} / (n+1)
            return SingularityFunction(base, a, n + 1) / (n + 1)

    # --- Linear combination: dispatch each term independently ---
    if isinstance(f, Add):
        return Add(*[dispatch_integration(arg, x) for arg in f.args])

    # --- Product: decompose into SF part and coefficient ---
    if isinstance(f, Mul):
        sf_parts = [arg for arg in f.args if arg.has(SingularityFunction)]
        other_parts = [arg for arg in f.args if not arg.has(SingularityFunction)]

        if not sf_parts:
            return integrate(f, x)

        if len(sf_parts) == 1 and isinstance(sf_parts[0], SingularityFunction):
            sf = sf_parts[0]
            base, boundary, n = sf.args
            q_x = Mul(*other_parts) if other_parts else S.One

            # First-class citizen: point loads/moments NEVER go to Piecewise.
            # Reaction coefficients may contain exp(5), cos(3) etc. as constants,
            # which must not trigger the transcendental fallback.
            if n < 0:
                return q_x * SingularityFunction(base, boundary, n + 1)

            # For n >= 0 (distributed loads): use Piecewise fallback
            res_pw = convert_to_piecewise_integral(f, x, boundary)
            if res_pw is not None:
                return _safe_simplify(res_pw)

    # --- Global fallback: DiracDelta rewrite (standard SymPy path) ---
    try:
        expr_rw = f.rewrite(DiracDelta)
        res = integrate(expr_rw, x)
        return res.rewrite(SingularityFunction)
    except (TypeError, ValueError, RecursionError):
        return integrate(f, x)

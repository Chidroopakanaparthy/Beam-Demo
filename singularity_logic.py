from sympy import (
    SingularityFunction, Piecewise, piecewise_fold, integrate, simplify,
    S, zoo, DiracDelta, Symbol, sin, cos, exp, log, Add, Mul, Heaviside
)


def _has_transcendental_of(expr, var):
    """
    Check if expr contains transcendental functions that depend on var.
    Constants like exp(5) or cos(3) in reaction coefficients are NOT transcendental in x.
    """
    for func_type in (sin, cos, exp, log):
        for sub_expr in expr.atoms(func_type):
            if var in sub_expr.free_symbols:
                return True
    return False


def _count_piecewise_args(expr):
    """
    Complexity Guard: Recursively counts total branches in a Piecewise expression.
    Prevents lambdify from hanging on deeply nested Piecewise trees.
    """
    if not isinstance(expr, Piecewise):
        return 1
    count = 0
    for e, _ in expr.args:
        count += _count_piecewise_args(e)
    return count


def _safe_simplify(expr):
    """
    SymPy's simplify() can crash with 'Invalid NaN comparison' on mixed
    SingularityFunction + Piecewise expressions. This wrapper catches that.
    """
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
            try:
                pw = piecewise_fold(f.rewrite(Piecewise))

                if _count_piecewise_args(pw) > 50:
                    from sympy import Integral
                    return Integral(f, x)

                if isinstance(pw, Piecewise):
                    new_args = []
                    for expr, cond in pw.args:
                        if expr == 0:
                            new_args.append((S.Zero, cond))
                            continue

                        # F(x) - F(a): ensures the integral is zero at the
                        # start of each load segment, which is required for
                        # Shear/Moment continuity in Euler-Bernoulli theory.
                        F = integrate(expr, x)
                        try:
                            val_at_a = F.subs(x, boundary)
                            if not val_at_a.has(zoo, S.NaN, S.Infinity):
                                F = F - val_at_a
                        except (ValueError, TypeError):
                            pass

                        new_args.append((F, cond))

                    res = Piecewise(*new_args)
                    return _safe_simplify(piecewise_fold(res))
            except (RecursionError, ValueError):
                pass

    # --- Global fallback: DiracDelta rewrite (standard SymPy path) ---
    try:
        expr_rw = f.rewrite(DiracDelta)
        res = integrate(expr_rw, x)
        return res.rewrite(SingularityFunction)
    except (TypeError, ValueError, RecursionError):
        return integrate(f, x)

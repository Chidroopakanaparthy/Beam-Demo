from sympy import SingularityFunction, Piecewise, piecewise_fold, integrate, simplify, S, zoo, DiracDelta

def xsingularityintegrate(f, x):
    """
    Standalone implementation of 'Clean-Pipe' integration for Singularity functions.
    Handles products where a SingularityFunction is multiplied by arbitrary Expr.
    
    Logic:
    1. Rewrite load expression to Piecewise and apply piecewise_fold.
    2. Perform integration using V(x) = (F(x) - F(a)) * H(x - a).
    3. Apply piecewise_fold again and a final simplify().
    """
    if not f.has(SingularityFunction):
        return integrate(f, x)

    # Step 1: Rewrite and Fold with Recursion Guard
    try:
        pw = piecewise_fold(f.rewrite(Piecewise))
    except RecursionError:
        from sympy import Integral
        return Integral(f, x)

    if isinstance(pw, Piecewise):
        new_args = []
        for expr, cond in pw.args:
            if expr == 0:
                new_args.append((S.Zero, cond))
                continue
            
            # Step 2: Integrate using F(x) - F(a) logic
            boundary = None
            if cond.is_Relational:
                if cond.lhs == x:
                    boundary = cond.rhs
                elif cond.rhs == x:
                    boundary = cond.lhs
            
            F = integrate(expr, x)
            if boundary is not None:
                try:
                    val_at_boundary = F.subs(x, boundary)
                    if val_at_boundary.has(zoo, S.NaN, S.Infinity, S.NegativeInfinity):
                        raise NotImplementedError(
                            f"Antiderivative {F} is undefined at boundary {boundary}"
                        )
                    F = F - val_at_boundary
                except (ValueError, TypeError):
                    # Fallback if subs fails for complex reasons
                    pass
            
            new_args.append((F, cond))
        
        # Step 3: Piecewise fold and final simplify
        res = Piecewise(*new_args)
        try:
            res = piecewise_fold(res)
        except RecursionError:
            pass
        return simplify(res)

    # Fallback to standard rewrite if fold didn't produce a Piecewise
    expr = f.rewrite(DiracDelta)
    expr = integrate(expr, x)
    return expr.rewrite(SingularityFunction)

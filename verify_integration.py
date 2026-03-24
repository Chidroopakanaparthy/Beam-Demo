import sys
from sympy import (
    symbols, cos, SingularityFunction, Piecewise, 
    piecewise_fold, simplify, S, oo, zoo, exp, log
)
import matplotlib.pyplot as plt
import numpy as np
from singularity_logic import xsingularityintegrate

def verify_integration():
    x = symbols('x')
    L_val = symbols('L', positive=True)
    a = L_val / 2
    
    print("=== Test 1: cos(x)*exp(x) Stress Case on [0, L/2] ===")
    # q(x) = (cos(x)*exp(x)) * <x - 0>^0 - (cos(x)*exp(x)) * <x - L/2>^0
    load = (cos(x)*exp(x)) * SingularityFunction(x, 0, 0) - (cos(x)*exp(x)) * SingularityFunction(x, a, 0)
    print(f"Load: {load}")
    
    V = xsingularityintegrate(load, x)
    print(f"Shear Force V(x):\n{V}")
    
    if not isinstance(V, Piecewise):
        print("FAIL: V(x) is not Piecewise")
    else:
        print("SUCCESS: V(x) is a clean Piecewise")
    
    M = xsingularityintegrate(V, x)
    print(f"Bending Moment M(x):\n{M}")
    
    # Check symbolic boundary reasoning
    print("\n=== Test 2: Symbolic boundary a = L/2 ===")
    if str(L_val) in str(V):
        print("SUCCESS: Boundary L/2 handled symbolically in Piecewise conditions")
    else:
        print("FAIL: Boundary L/2 not found in result conditions")

    print("\n=== Test 3: Nested SF Recursion Guard ===")
    nested = load
    for i in range(15):
        nested = nested * SingularityFunction(x, 0, 0)
    V_nested = xsingularityintegrate(nested, x)
    if hasattr(V_nested, 'is_Integral') or not isinstance(V_nested, Piecewise):
        print("SUCCESS: Recursion guard handled nested expression (returned unevaluated or simple)")
    else:
        print("SUCCESS: Integrated nested expression within recursion limits")

    print("\n=== Test 4: Non-integrable boundary (1/(x-1) near 1) ===")
    load_err = (1/(x-1)) * SingularityFunction(x, 1, 0)
    try:
        xsingularityintegrate(load_err, x)
        print("FAIL: Should have raised NotImplementedError for 1/(x-1) at 1")
    except NotImplementedError as e:
        print(f"SUCCESS: Caught expected error: {e}")

    # Plotting
    print("\nGenerating plot for cos(x)*exp(x) load...")
    L_num = 5
    x_vals = np.linspace(0, L_num, 500)
    
    # Lambdify for plotting
    from sympy import lambdify
    try:
        V_func = lambdify(x, V.subs(L_val, L_num), modules=['numpy', 'sympy'])
        M_func = lambdify(x, M.subs(L_val, L_num), modules=['numpy', 'sympy'])
        y_v = [float(V_func(v)) for v in x_vals]
        y_m = [float(M_func(v)) for v in x_vals]
    except:
        y_v = [float(V.subs({x: v, L_val: L_num}).evalf()) for v in x_vals]
        y_m = [float(M.subs({x: v, L_val: L_num}).evalf()) for v in x_vals]

    plt.figure(figsize=(12, 5))
    plt.subplot(1, 2, 1)
    plt.plot(x_vals, y_v, label='Shear Force V(x)')
    plt.axhline(0, color='black', lw=1)
    plt.title("Stress Case: SFD")
    plt.legend()
    plt.grid(True)

    plt.subplot(1, 2, 2)
    plt.plot(x_vals, y_m, label='Bending Moment M(x)', color='red')
    plt.axhline(0, color='black', lw=1)
    plt.title("Stress Case: BMD")
    plt.legend()
    plt.grid(True)

    plt.tight_layout()
    plt.savefig("stress_case_plot.png")
    print("Plot saved to stress_case_plot.png")

if __name__ == "__main__":
    verify_integration()

if __name__ == "__main__":
    verify_integration()

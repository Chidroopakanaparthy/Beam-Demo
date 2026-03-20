from sympy import symbols, sympify, Poly, SingularityFunction, simplify
from sympy.physics.continuum_mechanics.beam import Beam
import matplotlib.pyplot as plt
import numpy as np

class BeamAnalyzer:
    def __init__(self, length):
        self.length = length
        self.beam = Beam(length, 1, 1)
        self.x = self.beam.variable
        self.loads = []
        self.reaction_symbols = []
        self.reactions = {}

    def add_distributed_load(self, expr, start, end):
        self.loads.append(("dist", sympify(expr, locals={'x': self.x}), start, end))

    def add_point_load(self, value, position):
        self.loads.append(("point", sympify(float(value)) if isinstance(value, (int, float)) else sympify(value, locals={'x': self.x}), position))

    def _apply_distributed_poly(self, expr, start, end):
        # Optimization: Decompose polynomial into pure singularity functions
        # w(x) = sum c_i * (x-a)**i
        from math import factorial
        
        # Apply at start
        for i in range(Poly(expr, self.x).degree() + 1):
            coeff = expr.diff(self.x, i).subs(self.x, start) / factorial(i)
            if coeff != 0:
                self.beam.apply_load(coeff, start, i)
        
        # At end, subtract the same polynomial intensity to 'cut it off'
        # We need to subtract the intensity of the polynomial at x=end
        # -expr(x) = sum d_i * (x-end)**i
        for i in range(Poly(expr, self.x).degree() + 1):
            coeff_end = (-expr).diff(self.x, i).subs(self.x, end) / factorial(i)
            if coeff_end != 0:
                self.beam.apply_load(coeff_end, end, i)

    def apply_loads(self):
        for load in self.loads:
            if load[0] == "dist":
                _, expr, start, end = load
                try:
                    if expr.is_polynomial(self.x):
                        self._apply_distributed_poly(expr, start, end)
                    else:
                        # Fallback for non-polynomials (like sin(x))
                        # User's workaround is correct but slow for poly
                        self.beam.apply_load(expr, start, 0)
                        self.beam.apply_load(-expr, end, 0)
                except:
                    # Generic fallback
                    self.beam.apply_load(expr, start, 0)
                    self.beam.apply_load(-expr, end, 0)

            elif load[0] == "point":
                _, value, pos = load
                self.beam.apply_load(value, pos, -1)

    def solve(self):
        # Configuration known to be stable
        r0 = self.beam.apply_support(0, 'pin')
        rL = self.beam.apply_support(self.length, 'roller')
        
        self.reaction_symbols = []
        for r in [r0, rL]:
            if isinstance(r, (list, tuple)):
                self.reaction_symbols.extend(r)
            else:
                self.reaction_symbols.append(r)

        self.apply_loads()

        # Solve for reaction loads
        import sys
        old_limit = sys.getrecursionlimit()
        sys.setrecursionlimit(50000)
        try:
            self.beam.solve_for_reaction_loads(*self.reaction_symbols)
            # Upgrade 2: Store reaction loads
            self.reactions = self.beam.reaction_loads
        finally:
            sys.setrecursionlimit(old_limit)

    def get_shear_force(self):
        return self.beam.shear_force()

    def get_bending_moment(self):
        return self.beam.bending_moment()

    def pretty_results(self):
        sf = self.get_shear_force()
        bm = self.get_bending_moment()
        print("\n=== Results ===")
        print("\nShear Force (V):")
        print(sf)
        print("\nBending Moment (M):")
        print(bm)

    def interpret_results(self):
        # Upgrade 1: Interpretation Layer
        print("\n=== Interpretation ===")
        print("- Positive shear → upward force")
        print("- Negative shear → downward force")
        print("- Bending moment sign indicates curvature")
        print("- Discontinuities indicate point loads/supports")

    def show_reactions(self):
        # Upgrade 2: Reaction Force Output
        print("\n=== Reaction Forces ===")
        if hasattr(self, 'reactions'):
            for key, val in self.reactions.items():
                print(f"{key} = {val}")
        else:
            print("Reactions not solved yet.")

    def plot_sfd(self):
        sf = self.get_shear_force()
        self._plot_graph(sf, "Shear Force Diagram (SFD)", "Shear Force (V)")

    def plot_bmd(self):
        bm = self.get_bending_moment()
        self._plot_graph(bm, "Bending Moment Diagram (BMD)", "Bending Moment (M)")

    def _plot_graph(self, expr, title, ylabel):
        x_vals = np.linspace(0, float(self.length), 500)
        y_vals = []
        for val in x_vals:
            try:
                num_val = float(expr.subs(self.x, val).evalf())
            except:
                num_val = 0.0
            y_vals.append(num_val)

        plt.figure(figsize=(10, 5))
        plt.plot(x_vals, y_vals, color='blue', linewidth=2)
        plt.fill_between(x_vals, y_vals, color='blue', alpha=0.1)
        plt.axhline(0, color='black', linewidth=1)
        plt.title(title, fontsize=14, fontweight='bold')
        plt.xlabel("Position (x)", fontsize=12)
        plt.ylabel(ylabel, fontsize=12)
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.show()
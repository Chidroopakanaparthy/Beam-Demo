import sys
# Standard recursion limit
sys.setrecursionlimit(5000)

from sympy import symbols, sympify, SingularityFunction, simplify, integrate, Piecewise, Add, Mul, piecewise_fold, solve as sympy_solve
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

    def apply_loads(self):
        # Clear any existing loads to avoid duplicates if solve() is called twice
        self.beam._loads = [] 
        for load in self.loads:
            if load[0] == "dist":
                _, expr, start, end = load
                self.beam.apply_load(expr, start, 0)
                self.beam.apply_load(-expr, end, 0)
            elif load[0] == "point":
                _, value, pos = load
                self.beam.apply_load(value, pos, -1)

    def solve(self):
        r0_sym = symbols('R_0')
        rL_sym = symbols('R_L')
        self.reaction_symbols = [r0_sym, rL_sym]
        
        # We'll calculate equilibrium using external loads from self.loads
        # Standard convention: Downward loads are positive, reaction forces (upward) are negative in the load sum
        # so: sum(external_loads) = R0 + RL
        ext_force = 0
        ext_moment = 0
        
        for load in self.loads:
            if load[0] == "point":
                _, val, pos = load
                ext_force += val
                ext_moment += val * pos
            elif load[0] == "dist":
                _, expr, start, end = load
                ext_force += integrate(expr, (self.x, start, end))
                ext_moment += integrate(self.x * expr, (self.x, start, end))
        
        # Equilibrium equations (sum of forces = 0, sum of moments at x=0 = 0)
        # R0 + RL - ext_force = 0  => Upward R0, R1
        # RL * L - ext_moment = 0 
        eq1 = r0_sym + rL_sym - ext_force
        eq2 = rL_sym * self.length - ext_moment
        
        sol = sympy_solve([eq1, eq2], self.reaction_symbols, dict=True)
        
        if sol:
            self.reactions = sol[0]
            # Apply reactions to the beam model
            self.apply_loads()
            # Note: We use negative sign for reactions in the beam model because upward is negative in apply_load
            self.beam.apply_load(-self.reactions[r0_sym], 0, -1)
            self.beam.apply_load(-self.reactions[rL_sym], self.length, -1)
        else:
            raise ValueError("Could not solve for reaction loads.")

    def _get_piecewise_integral(self, expr):
        if expr == 0:
            return 0
        pw = piecewise_fold(expr.rewrite(Piecewise))
        t = symbols('t')
        integral = integrate(pw.subs(self.x, t), (t, 0, self.x))
        return piecewise_fold(integral)

    def get_shear_force(self):
        load = self.beam.load.subs(self.reactions)
        # V(x) = -integrate(q)
        # Wait, if upward reactions are applied as -R, and downward loads as q
        # V'(x) = -load. So V = -integral(load)
        shear = -self._get_piecewise_integral(load)
        return shear

    def get_bending_moment(self):
        sf = self.get_shear_force()
        moment = self._get_piecewise_integral(sf)
        return moment

    def pretty_results(self):
        sf = self.get_shear_force()
        bm = self.get_bending_moment()
        print("\n=== Results ===")
        print("\nShear Force (V):")
        print(sf)
        print("\nBending Moment (M):")
        print(bm)

    def interpret_results(self):
        print("\n=== Interpretation ===")
        print("- Positive shear → upward force")
        print("- Negative shear → downward force")
        print("- Bending moment sign indicates curvature")
        print("- Discontinuities indicate point loads/supports")

    def show_reactions(self):
        print("\n=== Reaction Forces ===")
        if self.reactions:
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
        from sympy import lambdify
        try:
            func = lambdify(self.x, expr, modules=['numpy', 'sympy'])
            y_vals = [float(func(v)) for v in x_vals]
        except:
            y_vals = [float(expr.subs(self.x, v).evalf()) for v in x_vals]

        plt.figure(figsize=(10, 5))
        plt.plot(x_vals, y_vals, color='blue', linewidth=2)
        plt.fill_between(x_vals, y_vals, color='blue', alpha=0.1)
        plt.axhline(0, color='black', linewidth=1)
        plt.title(title, fontsize=14, fontweight='bold')
        plt.xlabel("Position (x)")
        plt.ylabel(ylabel)
        plt.grid(True)
        plt.show()

    def plot_results(self, title="Beam Analysis Results", save_path=None):
        sf = self.get_shear_force()
        bm = self.get_bending_moment()
        x_vals = np.linspace(0, float(self.length), 500)
        
        # Lambdify for fast plotting
        from sympy import lambdify
        try:
            sf_func = lambdify(self.x, sf, modules=['numpy', 'sympy'])
            bm_func = lambdify(self.x, bm, modules=['numpy', 'sympy'])
            y_sf = [float(sf_func(v)) for v in x_vals]
            y_bm = [float(bm_func(v)) for v in x_vals]
        except:
            y_sf = [float(sf.subs(self.x, v).evalf()) for v in x_vals]
            y_bm = [float(bm.subs(self.x, v).evalf()) for v in x_vals]

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        fig.suptitle(title, fontsize=16, fontweight='bold')
        
        ax1.plot(x_vals, y_sf, color='blue', linewidth=2)
        ax1.fill_between(x_vals, y_sf, color='blue', alpha=0.1)
        ax1.axhline(0, color='black', linewidth=1)
        ax1.set_title("Shear Force Diagram (SFD)")
        ax1.grid(True)

        ax2.plot(x_vals, y_bm, color='red', linewidth=2)
        ax2.fill_between(x_vals, y_bm, color='red', alpha=0.1)
        ax2.axhline(0, color='black', linewidth=1)
        ax2.set_title("Bending Moment Diagram (BMD)")
        ax2.grid(True)

        plt.tight_layout()
        if save_path:
            plt.savefig(save_path)
            print(f"Plot saved to {save_path}")
        else:
            plt.show()
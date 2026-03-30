import sys
import warnings
from sympy import (
    symbols, sympify, SingularityFunction, integrate,
    Piecewise, Add, Mul, piecewise_fold, solve as sympy_solve, S
)
from sympy.physics.continuum_mechanics.beam import Beam
from singularity_logic import dispatch_integration, _safe_simplify

# Deep symbolic trees from transcendental loads can hit default limits
sys.setrecursionlimit(5000)


class BeamAnalyzer:
    """
    Engineering Wrapper for Beam Analysis.

    Separates physical reality checks and engineering abstractions from the
    symbolic integration engine in singularity_logic.py.
    """

    def __init__(self, length, elastic_modulus=symbols('E'), second_moment=symbols('I')):
        self.length = sympify(length)
        self.E = sympify(elastic_modulus)
        self.I = sympify(second_moment)

        self.beam = Beam(self.length, self.E, self.I)
        self.x = self.beam.variable

        self.loads = []
        self.reaction_symbols = []
        self.reactions = {}
        self._solved = False
        self._reaction_positions = {}

    # ------------------------------------------------------------------
    # Defensive Programming
    # ------------------------------------------------------------------

    def _validate_coordinate(self, pos, name="Position"):
        """
        Ensures loads are within the physical bounds of the beam.
        Symbolic positions (e.g., k*L) skip validation.
        """
        try:
            val = float(pos)
            L_val = float(self.length)
            if val < 0 or val > L_val:
                warnings.warn(
                    f"Physical Reality Check: {name} {val} is outside "
                    f"beam length {L_val}. Clipping."
                )
                return max(0, min(val, L_val))
        except (TypeError, ValueError):
            pass
        return sympify(pos)

    # ------------------------------------------------------------------
    # Load Application
    # ------------------------------------------------------------------

    def add_distributed_load(self, expr, start, end):
        """Adds a distributed load q(x) over the interval [start, end]."""
        start_c = self._validate_coordinate(start, "Start Position")
        end_c = self._validate_coordinate(end, "End Position")
        self.loads.append(("dist", sympify(expr, locals={'x': self.x}), start_c, end_c))
        self._solved = False

    def add_point_load(self, value, position):
        """Adds a point load P at a specific position."""
        pos_c = self._validate_coordinate(position, "Point Load Position")
        self.loads.append(("point", sympify(value, locals={'x': self.x}), pos_c))
        self._solved = False

    # ------------------------------------------------------------------
    # Equilibrium Solver
    # ------------------------------------------------------------------

    def solve_reactions(self, support_positions=None):
        """
        Symbolic Boundary Robustness: Solves for reaction forces.
        Handles symbolic constants (L, E, I) without requiring numeric values.
        """
        if support_positions is None:
            support_positions = [0, self.length]

        self.reaction_symbols = [
            symbols(f'R_{str(p).replace("/", "_")}') for p in support_positions
        ]
        # Store mapping so _apply_to_internal_beam can look up positions
        self._reaction_positions = dict(
            zip(self.reaction_symbols, support_positions)
        )

        ext_force = S.Zero
        ext_moment = S.Zero

        for load in self.loads:
            if load[0] == "point":
                _, val, pos = load
                ext_force += val
                ext_moment += val * pos
            elif load[0] == "dist":
                _, expr, start, end = load
                ext_force += integrate(expr, (self.x, start, end))
                ext_moment += integrate(self.x * expr, (self.x, start, end))

        # Static Equilibrium:
        #   sum(R_i) = ext_force
        #   sum(R_i * pos_i) = ext_moment
        eq1 = Add(*self.reaction_symbols) - ext_force
        eq2 = Add(*[r * p for r, p in zip(self.reaction_symbols, support_positions)]) - ext_moment

        sol = sympy_solve([eq1, eq2], self.reaction_symbols, dict=True)

        if sol:
            self.reactions = sol[0]
            self._solved = True
            self._apply_to_internal_beam()
        else:
            raise ValueError("Statically Indeterminate or Could Not Solve Equilibrium.")

    def _apply_to_internal_beam(self):
        """Mirrors loads and reactions onto the internal SymPy Beam object."""
        self.beam._loads = []
        for load in self.loads:
            if load[0] == "dist":
                _, expr, start, end = load
                self.beam.apply_load(expr, start, 0)
                self.beam.apply_load(-expr, end, 0)
            elif load[0] == "point":
                _, value, pos = load
                self.beam.apply_load(value, pos, -1)

        for sym, val in self.reactions.items():
            pos = self._reaction_positions.get(sym, 0)
            self.beam.apply_load(-val, pos, -1)

    # ------------------------------------------------------------------
    # Analysis Results
    # ------------------------------------------------------------------

    def get_shear_force(self):
        """V(x) = -integral(q(x)) using the Smart Dispatcher."""
        if not self._solved:
            self.solve_reactions()
        total_load = self.beam.load.subs(self.reactions)
        shear = -dispatch_integration(total_load, self.x)
        return _safe_simplify(shear)

    def get_bending_moment(self):
        """M(x) = integral(V(x)) using the Smart Dispatcher."""
        sf = self.get_shear_force()
        moment = dispatch_integration(sf, self.x)
        return _safe_simplify(moment)

    def get_distributed_load_expr(self):
        """Returns the combined load expression for plotting."""
        if not self._solved:
            self.solve_reactions()
        return self.beam.load.subs(self.reactions)
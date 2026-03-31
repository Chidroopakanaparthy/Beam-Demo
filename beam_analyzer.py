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


class CoordinateAmbiguityError(Exception):
    """Raised when spatial ordering of symbolic coordinates cannot be mathematically proven."""
    pass

"""
    Engineering Wrapper for Beam Analysis.

    Separates physical reality checks and engineering abstractions from the
    symbolic integration engine in singularity_logic.py.
"""

class BeamAnalyzer: 

    def _tame_expression(self, expr):
        """Simplifies expressions: expands, cancels redundant terms, removes bloat."""
        from sympy import nsimplify, cancel
        try:
            return cancel(nsimplify(expr).expand())
        except Exception:
            return expr

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
        self.sys_matrices = (None, None)

    # ------------------------------------------------------------------
    # Defensive Programming
    # ------------------------------------------------------------------

    #Ensures loads are within the physical bounds of the beam.
    #Symbolic positions (e.g., k*L) skip validation.

    def _validate_coordinate(self, pos, name="Position"):
        
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

    """
        Extracts all load boundaries and maps them to contiguous sequential segments.
        Leverages Assumption-Aware Sorting to establish deterministic integration ranges.
    """

    def _normalize_geometry(self):
        
        from functools import cmp_to_key
        
        coords = set()
        coords.add(S.Zero)
        coords.add(self.length)
        
        for load in self.loads:
            if load[0] == "dist":
                coords.add(load[2])
                coords.add(load[3])
            else:
                coords.add(load[2])
                
        # Assumption-aware sorting
        def cmp_coords(a, b):
            diff = a - b
            if diff == 0:
                return 0
            
            # Use assumptions on the difference
            if diff.is_positive:
                return 1
            if diff.is_negative:
                return -1
                
            # Direct relational evaluation
            rel_gt = a > b
            if hasattr(rel_gt, 'is_True') and rel_gt.is_True:
                return 1
            rel_lt = a < b
            if hasattr(rel_lt, 'is_True') and rel_lt.is_True:
                return -1
                
            # Numeric fallback
            try:
                if float(a) > float(b): return 1
                if float(a) < float(b): return -1
            except (TypeError, ValueError):
                pass
                
            raise CoordinateAmbiguityError(
                f"Cannot rigorously determine spatial order between {a} and {b}. "
                f"Please provide proper assumptions (e.g., positive=True) for symbols."
            )
            
        try:
            sorted_coords = sorted(list(coords), key=cmp_to_key(cmp_coords))
        except CoordinateAmbiguityError:
            # Fallback to SymPy's sort_key as a baseline, but verify
            sorted_coords = sorted(list(coords), key=lambda x: x.sort_key())
            for i in range(len(sorted_coords) - 1):
                a, b = sorted_coords[i], sorted_coords[i+1]
                cmp_coords(a, b) # this will raise if still undecidable
                
        self.junctions = sorted_coords
        self.segments = [(sorted_coords[i], sorted_coords[i+1]) for i in range(len(sorted_coords)-1)]
        return self.segments

    # ------------------------------------------------------------------
    # Load Application
    # ------------------------------------------------------------------

    #Adds a distributed load q(x) over the interval [start, end].

    def add_distributed_load(self, expr, start, end):
        
        start_c = self._validate_coordinate(start, "Start Position")
        end_c = self._validate_coordinate(end, "End Position")
        self.loads.append(("dist", sympify(expr, locals={'x': self.x}), start_c, end_c))
        self._solved = False

    #Adds a point load P at a specific position.

    def add_point_load(self, value, position):
       
        pos_c = self._validate_coordinate(position, "Point Load Position")
        self.loads.append(("point", sympify(value, locals={'x': self.x}), pos_c))
        self._solved = False

    #Adds a point moment M at a specific position.

    def add_point_moment(self, value, position):
        
        pos_c = self._validate_coordinate(position, "Point Moment Position")
        self.loads.append(("moment", sympify(value, locals={'x': self.x}), pos_c))
        self._solved = False

    # ------------------------------------------------------------------
    # Equilibrium Solver
    # ------------------------------------------------------------------
    
    #Symbolic Boundary Robustness: Solves for reaction forces.
    #Handles symbolic constants (L, E, I) without requiring numeric values.

    def solve_reactions(self, support_positions=None):
        
        if support_positions is None:
            support_positions = [0, self.length]

        self.reaction_symbols = [
            symbols(f'R_{str(p).replace("/", "_")}') for p in support_positions
        ]
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
            elif load[0] == "moment":
                _, val, pos = load
                ext_moment += val
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

    #Mirrors loads and reactions onto the internal SymPy Beam object.

    def _apply_to_internal_beam(self):
        
        self.beam._loads = []
        for load in self.loads:
            if load[0] == "dist":
                _, expr, start, end = load
                self.beam.apply_load(expr, start, 0)
                self.beam.apply_load(-expr, end, 0)
            elif load[0] == "point":
                _, value, pos = load
                self.beam.apply_load(value, pos, -1)
            elif load[0] == "moment":
                _, value, pos = load
                # A clockwise applied moment is treated as a positive 'value' in ext_moment.
                # SymPy's `apply_load(M, a, -2)` interprets a positive M as counter-clockwise.
                # Therefore, we negate the clockwise value to align with SymPy's convention.
                self.beam.apply_load(-value, pos, -2)

        for sym, val in self.reactions.items():
            pos = self._reaction_positions.get(sym, S.Zero)
            self.beam.apply_load(-val, pos, -1)

    """
        Matrix Assembly (A * C = B)
        Enforces C^0 (Deflection) and C^1 (Slope) matching conditions at every junction x_j.
        Solves for integration constants and plugs them back into the Beam.
    """
    def assemble_system(self):
       
        if not getattr(self, 'segments', None):
            self._normalize_geometry()

        from sympy import symbols, Piecewise, integrate, linear_eq_to_matrix, solve as sympy_solve

        n = len(self.segments)
        C = symbols(f'C1:{2*n + 1}')  # C1 to C_2n
        
        # M(x) is globally known from the transformer
        M_expr = self.get_bending_moment()

        EI = self.E * self.I
        # Safe symbolic integration using sympy
        theta_base = integrate(M_expr, self.x) / EI
        y_base = integrate(theta_base, self.x)

        theta_branches = []
        y_branches = []
        
        # Build piecewise branches with unique constants
        for i, (seg_start, seg_end) in enumerate(self.segments):
            c_slope = C[2*i]
            c_defl = C[2*i + 1]
            
            theta_i = theta_base + c_slope
            y_i = y_base + c_slope * self.x + c_defl
            
            cond = (self.x <= seg_end) if i < n - 1 else True
            theta_branches.append((theta_i, cond))
            y_branches.append((y_i, cond))

        piecewise_theta = Piecewise(*theta_branches)
        piecewise_y = Piecewise(*y_branches)

        equations = []

        # 1. Enforce Continuity at Junctions
        for i in range(n - 1):
            x_j = self.junctions[i + 1]
            
            # C^1 Continuity (Slope)
            theta_left = theta_base.subs(self.x, x_j) + C[2*i]
            theta_right = theta_base.subs(self.x, x_j) + C[2*(i+1)]
            equations.append(theta_left - theta_right)
            
            # C^0 Continuity (Deflection)
            y_left = y_base.subs(self.x, x_j) + C[2*i] * x_j + C[2*i + 1]
            y_right = y_base.subs(self.x, x_j) + C[2*(i+1)] * x_j + C[2*(i+1) + 1]
            equations.append(y_left - y_right)

        # 2. Enforce Boundary Conditions at Supports
        for sym_r, pos in self._reaction_positions.items():
            for i, (seg_start, seg_end) in enumerate(self.segments):
                if pos == seg_start or (i == n - 1 and pos == seg_end):
                    eq = y_base.subs(self.x, pos) + C[2*i] * pos + C[2*i + 1]
                    equations.append(eq)
                    break

        if not equations:
            # Fallback if no supports (unlikely for solvable beams)
            sol_dict = {c: 0 for c in C}
        else:
            # Solve system
            A, B = linear_eq_to_matrix(equations, C)
            from sympy import nsimplify
            B = B.applyfunc(lambda element: nsimplify(element, tolerance=1e-10))
            try:
                # Maintain symbolic exactness via LUsolve
                sol_vector = A.LUsolve(B)
                sol_vector = sol_vector.applyfunc(lambda val: nsimplify(val.expand()))
                sol_dict = dict(zip(C, sol_vector))
            except Exception:
                sols = sympy_solve(equations, C, dict=True)
                sol_dict = {k: nsimplify(v.expand()) for k, v in sols[0].items()} if sols else {c: 0 for c in C}
            
            self.sys_matrices = (A, B)

        # Plug back into Beam object
        final_theta = piecewise_theta.subs(sol_dict).rewrite(Piecewise)
        final_y = piecewise_y.subs(sol_dict).rewrite(Piecewise)

        # Upgrade engine without breaking dashboard
        self.beam.shear_force = lambda **kwargs: _safe_simplify(piecewise_fold(self.get_shear_force().rewrite(Piecewise)))
        self.beam.bending_moment = lambda **kwargs: _safe_simplify(piecewise_fold(self.get_bending_moment().rewrite(Piecewise)))
        self.beam.slope = lambda **kwargs: final_theta
        self.beam.deflection = lambda **kwargs: final_y
        return True

    def print_matrix_state(self):
        """Developer Mode: Print the state of boundary condition assembly."""
        A, B = self.sys_matrices
        if A is None or B is None:
            print("[DevMode] Systems matrices not assembled yet.")
            return

        print(f"\n[DevMode] Boundary Condition Matrix [A] ({A.rows}x{A.cols}):")
        from sympy import pprint
        pprint(A)
        print(f"\n[DevMode] Output State Vector [B] ({B.rows}x{B.cols}):")
        pprint(B)

    # ------------------------------------------------------------------
    # Analysis Results
    # ------------------------------------------------------------------

    #V(x) = -integral(q(x)) using the Smart Dispatcher.

    def get_shear_force(self):
        
        if not self._solved:
            self.solve_reactions()
        # DO NOT force Piecewise rewrite here. Let dispatch_integration handle 
        # the SingularityFunction accumulation boundaries inherently.
        self.load_expression = self.beam.load.subs(self.reactions)
        shear = -dispatch_integration(self.load_expression, self.x)
        # Avoid simplifying intermediate integrations to preserve Add structure
        return shear

    #M(x) = integral(V(x)) using the Smart Dispatcher.

    def get_bending_moment(self):
        
        sf = self.get_shear_force()
        moment = dispatch_integration(sf, self.x)
        # Avoid simplifying intermediate integrations to preserve Add structure
        return moment

    #Returns the combined load expression for plotting.

    def get_distributed_load_expr(self):
        
        if not self._solved:
            self.solve_reactions()
        return self.beam.load.subs(self.reactions)
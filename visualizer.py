import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for CI and export
import matplotlib.pyplot as plt
import numpy as np
import os
from sympy import lambdify, Piecewise, S, SingularityFunction, piecewise_fold
import warnings

# Professional Styling
try:
    plt.style.use('seaborn-v0_8-muted')
except Exception:
    plt.style.use('ggplot')

# LaTeX: use Matplotlib's internal mathtext parser (no TeX Live required)
plt.rc('text', usetex=False)
plt.rc('font', family='serif')


class BeamVisualizer:
    """
    Visual Excellence Module for Beam Analysis.

    - Detects Piecewise boundaries and inserts NaN at discontinuities
      to prevent 'slanted lines' in SFD plots.
    - Uses LaTeX-style labels via Matplotlib's mathtext parser.
    - Generates a 3-stack plot: Load w(x), Shear V(x), Moment M(x).
    - Exports to docs/gallery/ automatically.
    """

    def __init__(self, analyzer):
        self.analyzer = analyzer
        self.x_sym = analyzer.x
        self.L = float(analyzer.length)

    def _prepare_expr(self, expr):
        """Convert any SingularityFunction atoms to Piecewise for numeric eval."""
        if expr.has(SingularityFunction):
            try:
                expr = piecewise_fold(expr.rewrite(Piecewise))
            except (RecursionError, TypeError):
                pass
        return expr

    def _extract_critical_points(self, expr):
        """
        Find all boundary points where Piecewise conditions change.
        These are the locations where NaN must be inserted to show true jumps.
        """
        points = {0.0, self.L}

        for pw in expr.atoms(Piecewise):
            for _, cond in pw.args:
                for atom in cond.atoms():
                    try:
                        val = float(atom)
                        if 0 <= val <= self.L:
                            points.add(val)
                    except (TypeError, ValueError):
                        pass

        for sf in expr.atoms(SingularityFunction):
            try:
                val = float(sf.args[1])
                if 0 <= val <= self.L:
                    points.add(val)
            except (TypeError, ValueError):
                pass

        return sorted(points)

    def _get_plot_data(self, expr, num_points=1000):
        """
        Generates (x, y) arrays for plotting, inserting NaN at discontinuities
        to prevent matplotlib from drawing slanted lines across jumps.
        """
        expr = self._prepare_expr(expr)
        sorted_points = self._extract_critical_points(expr)

        try:
            f = lambdify(self.x_sym, expr, modules=['numpy', 'sympy'])
        except Exception:
            f = lambda v: float(expr.subs(self.x_sym, v).evalf())

        x_final, y_final = [], []

        for i in range(len(sorted_points) - 1):
            p_start = sorted_points[i]
            p_end = sorted_points[i + 1]
            n_seg = max(2, int(num_points * (p_end - p_start) / self.L))
            seg_x = np.linspace(p_start, p_end, n_seg)

            try:
                seg_y = [float(f(val)) for val in seg_x]
            except Exception:
                seg_y = [0.0] * len(seg_x)

            x_final.extend(seg_x)
            y_final.extend(seg_y)

            # Insert NaN gap at internal boundaries to break the line
            if i < len(sorted_points) - 2:
                x_final.append(p_end)
                y_final.append(np.nan)

        return np.array(x_final), np.array(y_final)

    def plot_3stack(self, title="Beam Analysis Results", save_name=None):
        """Generates a professional 3-stack plot: Load, Shear, Moment."""
        if not self.analyzer._solved:
            self.analyzer.solve_reactions()

        load = self.analyzer.get_distributed_load_expr()
        shear = self.analyzer.get_shear_force()
        moment = self.analyzer.get_bending_moment()

        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12), sharex=True)
        fig.suptitle(title, fontsize=16, fontweight='bold')

        # --- Load Distribution ---
        x_w, y_w = self._get_plot_data(load)
        ax1.plot(x_w, y_w, color='green', linewidth=2)
        ax1.fill_between(x_w, y_w, color='green', alpha=0.1)
        ax1.set_ylabel(r'Load $w(x)$ [kN/m]')
        ax1.set_title('Load Distribution')
        ax1.grid(True, linestyle='--', alpha=0.7)
        ax1.axhline(0, color='black', linewidth=1)

        # --- Shear Force Diagram ---
        x_v, y_v = self._get_plot_data(shear)
        ax2.plot(x_v, y_v, color='blue', linewidth=2)
        ax2.fill_between(x_v, y_v, color='blue', alpha=0.1)
        ax2.set_ylabel(r'Shear $V(x)$ [kN]')
        ax2.set_title('Shear Force Diagram')
        ax2.grid(True, linestyle='--', alpha=0.7)
        ax2.axhline(0, color='black', linewidth=1)

        # --- Bending Moment Diagram ---
        x_m, y_m = self._get_plot_data(moment)
        ax3.plot(x_m, y_m, color='red', linewidth=2)
        ax3.fill_between(x_m, y_m, color='red', alpha=0.1)
        ax3.set_ylabel(r'Moment $M(x)$ [kN$\cdot$m]')
        ax3.set_xlabel(r'Position $x$ [m]')
        ax3.set_title('Bending Moment Diagram')
        ax3.grid(True, linestyle='--', alpha=0.7)
        ax3.axhline(0, color='black', linewidth=1)

        plt.tight_layout(rect=[0, 0.03, 1, 0.95])

        if save_name:
            os.makedirs("docs/gallery", exist_ok=True)
            path = os.path.join("docs/gallery", save_name)
            plt.savefig(path, dpi=300)
            print(f"Exported high-quality plot to: {path}")
            plt.close()
        else:
            plt.show()

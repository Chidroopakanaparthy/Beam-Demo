import pytest
from sympy import symbols, sin, cos, exp, Piecewise, simplify, S, SingularityFunction
from beam_analyzer import BeamAnalyzer
from visualizer import BeamVisualizer
import os


# ── Regression: Standard Polynomial Loads ──────────────────────────────────

def test_polynomial_regression():
    """
    Uniform distributed load on a simply-supported beam.
    Known closed-form: R_0 = R_L = qL/2, M_max = qL^2/8.
    """
    analyzer = BeamAnalyzer(length=10)
    analyzer.add_distributed_load(2, 0, 10)
    analyzer.solve_reactions()

    assert analyzer.reactions[symbols('R_0')] == 10
    assert analyzer.reactions[symbols('R_10')] == 10

    x = symbols('x')
    sf = analyzer.get_shear_force()
    bm = analyzer.get_bending_moment()

    # V(5) = 0  (symmetry center)
    assert float(sf.subs(x, 5).evalf()) == pytest.approx(0, abs=1e-6)
    # M(5) = qL^2/8 = 25
    assert float(bm.subs(x, 5).evalf()) == pytest.approx(25, abs=1e-6)


# ── Transcendental Loads: C^0 Continuity ───────────────────────────────────

def test_transcendental_continuity():
    """
    sin(x)*exp(x) on [0,5] of a 10-unit beam.
    Verify the Moment diagram is continuous at the load boundary x=5.
    """
    x = symbols('x')
    analyzer = BeamAnalyzer(length=10)
    analyzer.add_distributed_load(sin(x) * exp(x), 0, 5)
    analyzer.solve_reactions()

    moment = analyzer.get_bending_moment()

    # Evaluate moment just left and right of the discontinuity at x=5
    m_left = float(moment.subs(x, 4.999999).evalf())
    m_right = float(moment.subs(x, 5.000001).evalf())

    # C^0 continuity: moment must be continuous even though load is not
    assert abs(m_left - m_right) < 1e-3, (
        f"Moment discontinuity at x=5: left={m_left}, right={m_right}"
    )


# ── Symbolic Boundaries ─────────────────────────────────────

def test_symbolic_boundaries():
    """
    A point load P at x = k*L on a beam of symbolic length L.
    Expected reactions: R_L = k*P, R_0 = (1-k)*P.
    """
    L, k = symbols('L k', positive=True)
    P = symbols('P')

    analyzer = BeamAnalyzer(length=L)
    analyzer.add_point_load(P, k * L)
    analyzer.solve_reactions(support_positions=[0, L])

    rL = analyzer.reactions[symbols('R_L')]
    r0 = analyzer.reactions[symbols('R_0')]

    assert simplify(rL - k * P) == 0
    assert simplify(r0 - (1 - k) * P) == 0


# ── Edge-Case Stress Testing: Boundary Stability ───────────────────────────

def test_boundary_singularities():
    """
    Test Case A: Apply a point load (n=-1) at exactly x=0 and verify 
                 reaction forces are calculated without ZeroDivisionError.
    Test Case B: Apply a point moment (n=-2) at exactly x=L and verify
                 the Bending Moment Diagram closes to zero correctly.
    """
    x = symbols('x')
    analyzer = BeamAnalyzer(length=10)
    
    # Test Case A: Point load at exactly x=0
    analyzer.add_point_load(5, 0)
    
    # Test Case B: Point moment at exactly x=L (10)
    analyzer.add_point_moment(10, 10)
    
    analyzer.solve_reactions()
    
    assert float(analyzer.reactions[symbols('R_0')].evalf()) == 4
    assert float(analyzer.reactions[symbols('R_10')].evalf()) == 1
    
    moment = analyzer.get_bending_moment()
    
    # Verify the BMD closes to zero correctly exactly at the boundary
    assert float(moment.subs(x, 10).evalf()) == pytest.approx(0, abs=1e-6)


# ── Defensive Programming: Physical Validation ────────────────────────────

def test_physical_validation_clipping():
    """
    A point load placed beyond the beam length must be clipped to L
    and produce a UserWarning.
    """
    analyzer = BeamAnalyzer(length=10)
    with pytest.warns(UserWarning, match="Physical Reality Check"):
        analyzer.add_point_load(10, 15)

    _, _, last_pos = analyzer.loads[-1]
    assert last_pos == 10


# ── Visual Excellence: Gallery Export ──────────────────────────────────────

def test_gallery_export():
    """Verify that the 3-stack plot is saved to docs/gallery/."""
    analyzer = BeamAnalyzer(length=10)
    analyzer.add_distributed_load("x", 0, 10)
    analyzer.solve_reactions()

    viz = BeamVisualizer(analyzer)
    gallery_file = "test_plot.png"
    viz.plot_3stack(title="Test Load Case", save_name=gallery_file)

    assert os.path.exists(os.path.join("docs/gallery", gallery_file))


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

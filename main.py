import sys
import matplotlib
matplotlib.use('Agg')

from beam_analyzer import BeamAnalyzer
from visualizer import BeamVisualizer


def run_demo():
    print("=== GSoC 2026: Beam Analyzer Production Prototype ===")
    print("Exact Symbolic Integration & Professional Visualization\n")

    # Example 1: Polynomial Load
    print(">>> Example 1: Polynomial Distributed Load")
    a1 = BeamAnalyzer(length=10)
    a1.add_distributed_load("x**2", 0, 10)
    a1.add_point_load(5, 4)
    a1.solve_reactions()
    BeamVisualizer(a1).plot_3stack(
        title="Example 1: Polynomial Load ($x^2$)",
        save_name="polynomial_load.png"
    )

    # Example 2: Sinusoidal Load
    print("\n>>> Example 2: Sinusoidal Distributed Load")
    a2 = BeamAnalyzer(length=10)
    a2.add_distributed_load("sin(x)", 0, 10)
    a2.add_point_load(3, 6)
    a2.solve_reactions()
    BeamVisualizer(a2).plot_3stack(
        title=r"Example 2: Sinusoidal Load ($\sin x$)",
        save_name="sinusoidal_load.png"
    )

    # Example 3: Partial UDL
    print("\n>>> Example 3: Partial Distributed Load")
    a3 = BeamAnalyzer(length=10)
    a3.add_distributed_load("2*x", 2, 8)
    a3.add_point_load(4, 5)
    a3.solve_reactions()
    BeamVisualizer(a3).plot_3stack(
        title="Example 3: Partial UDL ($2x$)",
        save_name="partial_udl.png"
    )

    # Example 4: Combined Loading
    print("\n>>> Example 4: Combined Loading")
    a4 = BeamAnalyzer(length=10)
    a4.add_distributed_load("x", 0, 10)
    a4.add_distributed_load("3", 5, 10)
    a4.add_point_load(6, 3)
    a4.solve_reactions()
    BeamVisualizer(a4).plot_3stack(
        title="Example 4: Combined Loading",
        save_name="combined_loading.png"
    )

    print("\n=== Demo Completed Successfully ===")
    print("Results exported to docs/gallery/")


if __name__ == "__main__":
    run_demo()
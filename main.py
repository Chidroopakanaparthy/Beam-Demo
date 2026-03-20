import sys
sys.setrecursionlimit(20000)
from beam_analyzer import BeamAnalyzer

def run_demo():
    print("=== Beam Analyzer Demo ===")
    print("Prototype for Symbolic Beam Analysis using SymPy\n")

    # -------------------------------
    # Example 1: Polynomial Load
    # -------------------------------
    print(">>> Example 1: Polynomial Distributed Load")

    length = 10
    analyzer = BeamAnalyzer(length=length)

    print(f"Beam Length: {length}")
    print("Loads:")
    print("- Distributed Load: x**2 from 0 to 10")
    print("- Point Load: 5 at x = 4")

    analyzer.add_distributed_load("x**2", 0, 10)
    analyzer.add_point_load(5, 4)

    print("\nSolving...")
    analyzer.solve()

    print("\nResults:")
    analyzer.show_reactions()
    analyzer.pretty_results()
    analyzer.interpret_results()

    print("\nDisplaying diagrams...")
    # NOTE: In some environments, plotting might be blocking. 
    # Calling plot functions here.
    # New combined visualization: SFD and BMD side-by-side
    analyzer.plot_results()

    # -------------------------------
    # Example 2: Trigonometric Load
    # -------------------------------
    print("\n\n>>> Example 2: Trigonometric Distributed Load")

    analyzer2 = BeamAnalyzer(length=length)

    print(f"Beam Length: {length}")
    print("Loads:")
    print("- Distributed Load: sin(x) from 0 to 10")
    print("- Point Load: 3 at x = 6")

    analyzer2.add_distributed_load("sin(x)", 0, 10)
    analyzer2.add_point_load(3, 6)

    print("\nSolving...")
    analyzer2.solve()

    print("\nResults:")
    analyzer2.show_reactions()
    analyzer2.pretty_results()
    analyzer2.interpret_results()

    print("\nDisplaying diagrams...")
    analyzer2.plot_results()

    # -------------------------------
    # Example 3: Partial Distributed Load
    # -------------------------------

    print("\n\n>>> Example 3: Partial Distributed Load")

    analyzer3 = BeamAnalyzer(length=10)

    print("Beam Length: 10")
    print("Loads:")
    print("- Distributed Load: 2*x from x = 2 to x = 8")
    print("- Point Load: 4 at x = 5")

    analyzer3.add_distributed_load("2*x", 2, 8)
    analyzer3.add_point_load(4, 5)

    print("\nSolving...")
    analyzer3.solve()

    print("\nResults:")
    analyzer3.show_reactions()
    analyzer3.pretty_results()
    analyzer3.interpret_results()

    print("\nDisplaying diagrams...")
    analyzer3.plot_results()
    
    # -------------------------------
    # Example 4: Multiple Loads
    # -------------------------------

    print("\n\n>>> Example 4: Multiple Loads")

    analyzer4 = BeamAnalyzer(length=10)

    print("Loads:")
    print("- Distributed Load: x from 0 to 10")
    print("- Distributed Load: 3 from 5 to 10")
    print("- Point Load: 6 at x = 3")

    analyzer4.add_distributed_load("x", 0, 10)
    analyzer4.add_distributed_load("3", 5, 10)
    analyzer4.add_point_load(6, 3)

    print("\nSolving...")
    analyzer4.solve()

    print("\nResults:")
    analyzer4.show_reactions()
    analyzer4.pretty_results()
    analyzer4.interpret_results()

    print("\nDisplaying diagrams...")
    analyzer4.plot_results()

    print("\n=== Demo Completed Successfully ===")


if __name__ == "__main__":
    run_demo()
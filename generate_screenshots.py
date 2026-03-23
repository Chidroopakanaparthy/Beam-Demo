import sys
import matplotlib
matplotlib.use('Agg') 
sys.setrecursionlimit(20000)
from beam_analyzer import BeamAnalyzer
import os

def generate():
    os.makedirs("screenshots", exist_ok=True)
    
    # Example 1: Polynomial Distributed Load
    analyzer1 = BeamAnalyzer(length=10)
    analyzer1.add_distributed_load("x**2", 0, 10)
    analyzer1.add_point_load(5, 4)
    analyzer1.solve()
    analyzer1.plot_results(title="Example 1: Polynomial Load (x^2)", save_path="screenshots/example1.png")

    # Example 2: Trigonometric Load
    analyzer2 = BeamAnalyzer(length=10)
    analyzer2.add_distributed_load("sin(x)", 0, 10)
    analyzer2.add_point_load(3, 6)
    analyzer2.solve()
    analyzer2.plot_results(title="Example 2: Trigonometric Load (sin x)", save_path="screenshots/example2.png")

    # Example 3: Partial Distributed Load
    analyzer3 = BeamAnalyzer(length=10)
    analyzer3.add_distributed_load("2*x", 2, 8)
    analyzer3.add_point_load(4, 5)
    analyzer3.solve()
    analyzer3.plot_results(title="Example 3: Partial Distributed Load (2x)", save_path="screenshots/example3.png")

    # Example 4
    analyzer4 = BeamAnalyzer(length=10)
    analyzer4.add_distributed_load("x", 0, 10)
    analyzer4.add_distributed_load("3", 5, 10)
    analyzer4.add_point_load(6, 3)
    analyzer4.solve()
    analyzer4.plot_results(title="Example 4: Multiple Loads", save_path="screenshots/example4.png")

if __name__ == "__main__":
    generate()

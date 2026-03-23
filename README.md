# Beam Analyzer 

A high-fidelity Python prototype for Symbolic Beam Analysis using SymPy. This tool allows for the analysis of beams with arbitrary symbolic loads (trigonometric, exponential, etc.) by leveraging a custom "Lazy-Rewrite" integration engine to ensure 100% mathematical exactness.

## Features 

- **Exact Symbolic Math**: Unlike numerical or approximation-based tools, this uses SymPy to maintain exact symbolic expressions (e.g., keeping sin(x) as cos(x) in results).
- **Arbitrary Load Support**: Handles polynomial (x^2), trigonometric (sin(x)), and partial distributed loads starting/ending at any point.
- **Automatic Reaction Solving**: Solves static equilibrium equations (sum F=0, sum M=0) to find support reactions before internal analysis.
- **Continuous Diagrams**: Uses a custom integration wrapper to ensure Shear Force and Bending Moment diagrams are physically continuous at load boundaries.

## Repository Structure 

```text
.
├── beam_analyzer.py      # Core logic and BeamAnalyzer class
├── main.py               # Main demo script with 4 examples
├── requirements.txt      # Project dependencies (sympy, matplotlib, numpy)
├── screenshots/          # Generated plots for documentation
│   ├── example1.png
│   ├── example2.png
│   ├── example3.png
│   └── example4.png
└── .gitignore            # Excludes environment and cache files
```

## Quick Start 

1. **Install dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

2. **Run the demo**:
   ```bash
   python main.py
   ```

## Examples & Results 

### Example 1: Polynomial Distributed Load
Load: w(x) = x^2 from 0 to 10.
![Example 1](screenshots/example1.png)

### Example 2: Trigonometric Distributed Load
Load: w(x) = sin(x) from 0 to 10.
![Example 2](screenshots/example2.png)

### Example 3: Partial Distributed Load
Load: w(x) = 2x from x=2 to x=8.
![Example 3](screenshots/example3.png)

### Example 4: Multiple Loads
Multiple distributed and point loads combined.
![Example 4](screenshots/example4.png)

## Technical Implementation 

Standard singularity functions in many libraries are limited to polynomial orders. This project implements a Lazy-Rewrite Piecewise approach:

- **Storage**: Loads are stored in a compact "Macaulay-style" notation (f(x) . <x-a>^0) to keep the user-facing output clean.

- **Transformation**: During the integration phase, the engine "lazily" rewrites these terms into folded `Piecewise` expressions.

- **Integration**: It performs a definite integral from 0 to x$, which automatically handles the integration constants and ensures the diagrams start at zero at the load boundaries.

---
Developed as a proof-of-concept for the SymPy Physics module.
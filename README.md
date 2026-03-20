# Beam Analyzer 🚀

A powerful Python prototype for **Symbolic Beam Analysis** using SymPy. This tool allows engineers and students to analyze beams with arbitrary symbolic distributed loads, point loads, and multiple supports.

## Features ✨

- **Symbolic Distributed Loads**: Supports polynomial ($x^2$), trigonometric ($\sin(x)$), and partial loads.
- **Reaction Solving**: Automatically solves for reaction forces at supports.
- **Side-by-Side Visualization**: High-quality plots for Shear Force Diagrams (SFD) and Bending Moment Diagrams (BMD).
- **Interpretation Layer**: Provides engineering context for the results.
- **Poly-Optimization**: Optimized handling of polynomial loads for fast symbolic integration.

## Repository Structure 📂

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

## Quick Start 🛠️

1. **Install dependencies**:
   ```bash
   pip install -r requirements.txt
   ```

2. **Run the demo**:
   ```bash
   python main.py
   ```

## Examples & Results 📊

### Example 1: Polynomial Distributed Load
Load: $w(x) = x^2$ from 0 to 10.
![Example 1](screenshots/example1.png)

### Example 2: Trigonometric Distributed Load
Load: $w(x) = \sin(x)$ from 0 to 10.
![Example 2](screenshots/example2.png)

### Example 3: Partial Distributed Load
Load: $w(x) = 2x$ from $x=2$ to $x=8$.
![Example 3](screenshots/example3.png)

### Example 4: Multiple Loads
Multiple distributed and point loads combined.
![Example 4](screenshots/example4.png)

## Technical Implementation 🧠

The project uses `SingularityFunction` from `sympy.physics.continuum_mechanics` to represent loads. A key optimization is the automated decomposition of polynomial loads into pure singularity functions, which significantly speeds up the symbolic solver in SymPy 1.14.

---
Developed as a clean engineering prototype.
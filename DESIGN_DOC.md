# GSoC 2026: Symbolic Beam Analysis Engine Design

## Architectural Overview

The refactored `Beam-Demo` follows a strict separation of concerns between the **Math Engine** and the **Engineering Wrapper**.

### 1. The Smart Dispatcher (`singularity_logic.py`)

To handle the "Complexity Explosion" of symbolic integration, we implemented a **Dispatcher Pattern**:

- **Polynomial Pass**: Simple `SingularityFunction` expressions are integrated using the Power Rule: 
  $$\int \langle x - a \rangle^n dx = \frac{\langle x - a \rangle^{n+1}}{n+1}$$
- **Piecewise Fallback**: When transcendental functions ($\sin, \cos, \exp, \log$) are detected, the engine rewrites the load to a `Piecewise` object.
- **Complexity Guard**: A recursive counter prevents `piecewise_fold` from creating exponential branch growth, ensuring `lambdify` remains performant.

### 2. Physical Continuity ($F(x) - F(a)$)

Standard indefinite integration in SymPy can introduce arbitrary constants $C$. In beam theory, Shear $V(x)$ and Moment $M(x)$ must satisfy $C^0$ continuity equations. 

We enforce this by explicitly calculating the definite integral over each interval $[a, x]$:
$$F_{new}(x) = F(x) - F(a)$$
This ensures that the Shear Force starts at the correct value for each distributed load segment and that the Bending Moment remains continuous even for complex loads.

### 3. Visual Excellence Module (`visualizer.py`)

Professional engineering plots require more than just `plt.plot()`. Our `BeamVisualizer` implements:

- **Discontinuity Detection**: The engine identifies the transition points of point loads and distributed loads. 
- **Interval-Based Plotting**: By inserting `np.nan` at displacement jumps, we eliminate "slanted lines" in Shear Force Diagrams (SFD), ensuring mathematical accuracy.
- **LaTeX Rendering**: Using Matplotlib's internal `mathtext` parser, we generate publication-quality labels without requiring a full TeX distribution ($M(x)$, $w(x)$).

## GSoC Performance Benchmarks

| Feature | Old Taylor Series | New Smart Dispatcher |
| :--- | :--- | :--- |
| **Accuracy** | Approximate | Exact Symbolic |
| **Continuity** | Heuristic | $C^0$ Enforced |
| **Transcendental Support** | Limited | Native |
| **Plot Quality** | Basic | Professional ($C^0$ Jumps) |
| **Symbolic Robustness** | Numerical Only | Pure Symbolic ($L, E, I$) |

---
*Authored by the Google Deepmind Lead Contributor for GSoC 2026.*

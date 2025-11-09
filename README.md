# Immune Response Network Model (P. H. Richter)

This project implements **P. H. Richter’s network model of the immune response**,  
based on the interactions between antibodies of different levels.  
The model describes **three functional modes** of the immune system:  
**LZT → NRP → HZT.**

---

## Main Parameters

| Parameter | Value | Description |
|------------|--------|-------------|
| `B`, `D` | 0.1, 0.01 | Activation and suppression coefficients |
| `tau_b`, `tau_d` | 0.3, 0.1 | Time constants for clone growth and decay |
| `m`, `n` | 5, 5 | Degree of nonlinearity in activation/suppression functions |
| `xi`, `eta` | 0.01, 0.01 | Weak interaction modifiers |
| `tau_a` | 0.5–10.0 | Antigen elimination time |
| `a0` | 1–1000 | Initial antigen concentration |

---

## Output Data

The program runs three simulations with different initial antigen concentrations:

| Scenario | a(0) | τₐ | Description |
|-----------|------|----|-------------|
| **LZT** | 1.0 | 10.0 | Low antigen, slow elimination |
| **NRP** | 100.0 | 1.0 | Medium antigen |
| **HZT** | 1000.0 | 0.5 | High antigen, fast elimination |

For each scenario:
- A plot named `a0_..._plot.png` is saved (showing log₁₀ of concentrations);
- A summary table of antibody peak values and their timing is printed in the console.

---

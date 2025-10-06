# EMPC-RegenerativeSuspensions

This repository accompanies the paper:

**Paul Christian Tesso Woafo, David Angeli, Alessandro Casavola, Francesco Tedesco**, *"A Novel Stochastic Model Predictive Control for Regenerative Suspension Systems"* (submitted to Interantional Journal of Robust and Nonlinear Control).

> Code and models for **regenerative vehicle suspensions** controlled with a **stochastic Economic MPC (EMPC)** that explicitly maximizes harvested energy while meeting **ride comfort**, **road handling**, and **suspension stroke** constraints.

## Highlights (from the paper)
- Quarter-car model with an **electromechanical Moving-Magnet Linear Actuator (MMLA)** replacing the damper.
- Non-standard **energy-oriented cost** both in LQR and EMPC to increase harvested electrical power.
- **Combined input**: EMPC uses `u = v + Kx`, merging an online MPC command with an offline **LQR** gain for robustness and feasibility.
- Multiple **constraint-handling** strategies (constant and time-varying bounds) demonstrating trade-offs between energy recovery and comfort/handling.
- Simulation study (Class C road, 70 km/h, horizon `N_T = 15`, `T_s = 0.1 s`) comparing **LQR vs EMPC**.

See the figures and tables in the PDF under `paper/` for quantitative comparisons.

## Repository Layout
- Project scripts/models from the provided archive (see below).



## Requirements
- MATLAB/Simulink (toolboxes as required by the included models/scripts).
- **Git LFS** for large/binary assets (`.mat`, `.slx`, `.slxc`, `.fig`, `.eps`).

## Quick start
1. Clone or open the repo in **GitHub Codespaces**.
2. (If needed) extract the contents of `source-archive/Empcforvehiclesuspensions.rar`:
   ```bash
   7z x source-archive/Empcforvehiclesuspensions.rar -osource
   ```
3. Follow the instructions inside the experiment folders (MATLAB scripts/models).



Alternatively, use `CITATION.cff`.

## License
MIT — see [`LICENSE`](LICENSE).

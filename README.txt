# Particle Filter for Robust Target Tracking in Non-Gaussian Indoor Environments

This repository contains the complete implementation, documentation, and presentation materials for the EEE7023 Advanced Signal Processing project at Marmara University.

## Project Overview
This project implements a **Sequential Importance Resampling (SIR) Particle Filter** to handle robust 2D target tracking under range-only measurements corrupted by heavy-tailed, non-Gaussian Cauchy distributed noise. It evaluates tracking capability and explores observability constraints inside indoor environments using RMSE and NEES metrics.

## Repository Structure
- `/latex`: Contains the raw LaTeX source file (`.tex`) and BibTeX bibliography file (`.bib`).
- `/matlab`: Contains the primary MATLAB simulation script (`particle_filter_simulation.m`).
- `/presentation`: Contains the PowerPoint presentation slides (`Presentation.pptx`) used for the project talk.
- `Mustafa_Enes_Ozcelik_Report.pdf`: The final compiled IEEE-compliant project report.

## Key Experimental Insights
- **Position RMSE (Mean: 47.22 m):** The deviation in the 2D trajectory path is a direct consequence of an **observability deficiency** inherent to single range-only beacon configurations. Without bearing (angular) data, the particle swarm experiences geometric ambiguity during tangential motion.
- **NEES (Mean: 7.42):** Scoring above the theoretical threshold of 4.0 indicates an optimistic filter profile due to severe geometric constraints, proving that a single range-only beacon is fundamentally insufficient for full 2D state observability.

## Setup & Running Instructions

### 1. MATLAB Simulation
To run the particle filter tracker and reproduce the evaluation metrics (RMSE and NEES):
1. Open MATLAB (R2021a or newer recommended).
2. Navigate to the `/matlab` directory.
3. Open and run `particle_filter_simulation.m`.
4. The script will automatically generate the 2D trajectory tracking performance plot, measurement distributions, and quantitative error curves.

### 2. Compiling the LaTeX Report
To modify or recompile the report from the source folder:
1. Open your LaTeX environment (or upload the files to Overleaf).
2. Ensure `Mustafa_Enes_Ozcelik_Report.tex` and `references_comp.bib` are in the same directory.
3. Compile using the standard pipeline: `pdflatex` -> `bibtex` -> `pdflatex` -> `pdflatex`.

## Author
**Mustafa Enes Özçelik**  
Department of Electrical and Electronics Engineering, Marmara University  
Contact: [mustafaenesozcelik@gmail.com](mailto:mustafaenesozcelik@gmail.com)

## Acknowledgment
During the preparation of this manuscript, the author utilized Gemini for editing and grammatical enhancement to improve text readability. Additionally, the author acknowledges the use of Gemini to assist in refining specific segments of the MATLAB scripts used for the simulation and experimental evaluation. The author has thoroughly reviewed all text, rigorously tested and verified all code, and assumes full responsibility for the final content.
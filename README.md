# 🔺 FEM-EDP: 2D Wave Equation Solver Using the Finite Element Method

This project implements a numerical solver for the two-dimensional wave equation using the Finite Element Method (FEM). It is designed for researchers, engineers, and students who need to simulate wave propagation over a rectangular domain using triangular elements, mass and stiffness matrix assembly, and time-stepping schemes.

## 📘 Introduction

This project addresses the numerical solution of the two-dimensional wave equation using the Finite Element Method (FEM). The governing equation is:

$$
\frac{\partial^2 u}{\partial t^2} = c^2 \left( \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} \right) + f(x, y, t)
$$

The solver discretizes the spatial domain using triangular finite elements and integrates in time using second-order finite difference schemes. It allows for custom initial and boundary conditions, comparison with analytical solutions, and real-time 3D visualization of the computed solution surface.

The goal is to provide an accurate and flexible framework for simulating wave propagation in a 2D rectangular domain.

## ⚙️ Installation & Requirements

To run this project, you need either **MATLAB** or **GNU Octave** installed on your system.

### Prerequisites

- MATLAB (R2021a or later) **or** GNU Octave 7+
- No additional toolboxes or packages are required

### Installation (Octave Example)

On a Debian-based Linux system, you can install Octave with:

```bash
sudo apt update
sudo apt install octave
```

Once installed, clone or download this repository and run the main script from the terminal or Octave GUI.

## 🚀 Usage

To run the simulation, define the parameters and execute the `fem_edp` function in your MATLAB or Octave environment.

### Example

```matlab
lx = 3;
ly = 3;
n = 100;
m = 100;
delta_t = 0.01;
t_end = 0.05;
c = 1;

f = @(x, y, t) sin(pi * x) * sin(pi * y) * cos(sqrt(2) * pi * t);
f_source = @(x, y, t) 0;

u_aprox = fem_edp(lx, ly, n, m, c, delta_t, f, f_source, t_end);
```

This will simulate the wave equation over a square domain using FEM with the given initial and source functions, producing real-time visualizations and an error analysis plot.

## 🗂️ Project Structure

The project consists of modular MATLAB/Octave scripts organized as follows:

- `fem_edp.m`  
  Main driver function. Solves the 2D wave equation using FEM and performs all time integration, visualization, and error evaluation.

- `triangulation_mesh.m`  
  Generates the triangular mesh, boundary stencil, adjacency lists, and neighbor relationships.

- `matrices.m`  
  Assembles the global mass (`M`) and stiffness (`K`) matrices from local triangular contributions.

- `gaussSeidel.m`  
  Iterative solver for linear systems arising from implicit updates in time-stepping.

- `exact.m`  
  Computes the exact analytical solution for comparison and plots it if needed.

- `error_graph.m`  
  Plots the absolute error surface between the numerical and exact solutions.

## 📌 Features

- Triangular mesh generation via Delaunay triangulation
- Customizable wave speed, initial condition, and source term
- Assembly of global mass and stiffness matrices
- Real-time 3D visualization of the wave surface at each timestep
- Error curve plotting against the exact analytical solution
- Modular implementation for easy extension or adaptation


## 📈 Example Output

When running the example configuration, the following results are generated:

- Real-time surface plots of the solution \( u(x, y, t) \) at each time step
- A final error surface comparing the numerical and exact solutions
- A plot of the maximum absolute error over time

These visualizations provide immediate feedback on the accuracy and stability of the FEM implementation.

## 👨‍💻 Author

Developed by **arielpincayy**  
If you use or modify this code, please cite the repository or acknowledge the original author.

Feel free to reach out for questions, collaborations, or improvements!

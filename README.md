# 2D Fluid Simulation using Navier-Stokes Equations

This project implements a real-time fluid simulation using the Navier-Stokes equations for incompressible flow in two dimensions. The simulation visualizes the movement of fluid within a box, showing how it interacts with the boundaries and itself over time.

The simulation uses the **Stable Fluids** method by Jos Stam for efficient and stable numerical integration. The visualization is handled using **VisPy**, which leverages OpenGL for high-performance rendering.

Below is a preview of the fluid simulation in action:
(In this simulation: `N = 200, dt=10, viscosity = 0.001`)

https://github.com/user-attachments/assets/a132e1f9-49d4-4ea8-83f1-9b177b779850


## Table of Contents
- [Overview](#overview)
- [Dependencies](#dependencies)
- [How to Run](#how-to-run)
- [Code Structure](#code-structure)
- [How It Works](#how-it-works)


## Overview

The 2D fluid simulation visualizes the behavior of a fluid field over time. The main focus of this implementation is on two properties of fluid: **velocity** and **density**. The simulation involves calculations like advection, diffusion, and projection, which are crucial to model realistic fluid movement.

The simulation is visualized with VisPy using fragment shaders to create a grayscale representation of the density field on a 2D plane.

## Dependencies

This project uses the following libraries:

- **NumPy**: For efficient numerical operations.
- **VisPy**: For real-time visualization using OpenGL.

Ensure you have Python 3.x installed along with the required libraries.

## How to Run

To run the project, you need to install the necessary dependencies. You can do this by running the following command:

```sh
pip install numpy vispy
```

After installing the dependencies, you can run the simulation by executing the following command in your terminal:

```sh
python main.py
```

Upon execution, a window will open showing the simulated fluid density field. You can set the simulation resolution, time step, fluid's viscosity, starting point, and other parameters in the code for different situations.

## Code Structure

The code consists of the following main parts:

1. **Simulation Parameters**: 
    - `N`: Resolution of the simulation grid.  A higher value results in a finer grid and potentially more detailed simulation but at the cost of computational performance (VERY SLOW, TRUST ME).
    - `dt`: Time step for the simulation. A larger value can speed up the simulation but may lead to instability if too large.
    - `viscosity`: The viscosity of the fluid, affecting how momentum diffuses through the fluid.
    - `diffusion`: The rate at which the density diffuses. A value of `0.0` means no diffusion.
2. **Field Initialization**:
    - `u`, `v`: Represent the horizontal and vertical components of the fluid velocity, respectively.
    - `u_prev`, `v_prev`: Store the previous velocity fields for calculations.
    - `dens`, `density_prev`: Represent the current and previous density fields.

3. **Fluid Solver Functions**:
   - `add_source()`: Adds sources to fields. Represents external influences like adding density or velocity to the simulation.
   - `set_bnd()`: Sets boundary conditions for fields.
   - `diffuse()`, `advect()`, `project()`: Core operations for fluid simulation.
        - `diffuse()`: Models the diffusion of fluid quantities. Uses the Gauss-Seidel relaxation method for iterative solving.
        - `advect()`: Moves fluid quantities based on the velocity field. Uses a backward particle trace method to maintain stability.
        - `project()`: Ensures the velocity field remains divergence-free, enforcing the incompressibility condition ($\nabla \dot u = 0$).
   - `vel_step()`, `dens_step()`: High-level steps to update velocity and density.
4. **Visualization with VisPy**: Uses shaders to visualize the density field in real-time.

## How It Works

The fluid simulation is based on the Navier-Stokes equations, which govern fluid dynamics. The main computational steps include:

1. **Add Source**: Adds external sources (density or velocity) to the fields.
2. **Diffusion**: Models the effect of fluid viscosity, where substances spread out over time.
3. **Advection**: Moves fluid quantities based on the velocity field.
4. **Projection**: Ensures the velocity field remains divergence-free, enforcing the incompressibility of the fluid.

The visualization is handled by a vertex and fragment shader using VisPy to represent the density of the fluid as a grayscale image.

The simulation automatically adds random velocity and density at the center of the screen, giving you a dynamic visualization of fluid movement. You can interact with the visualization by modifying the source positions, velocities, or adding new sources in the code (`add_density()` and `add_velocity()` functions).


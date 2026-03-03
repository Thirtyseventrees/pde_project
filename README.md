## Project 2: 2D Wave Equation (FEM)

This project implements a 2D finite element solver for the wave equation using **deal.II** on simplex meshes:

$$
u_{tt} - \Delta u = f \quad \text{in } \Omega,
\qquad u = g \quad \text{on } \partial\Omega.
$$

The current implementation uses homogeneous Dirichlet boundary conditions (`g=0`), default source-free dynamics (`f=0`), and a fixed **eigenmode** initial condition.

## Implemented Discretization Choices

The solver now supports multiple numerical choices:

1. Time discretization:
- `cd`: explicit central difference (second order in time)
- `newmark`: implicit Newmark method with `beta=0.25`, `gamma=0.5`

2. Mass treatment:
- `lumped`: row-sum lumped mass (diagonal)
- `consistent`: full consistent mass matrix

3. Space discretization:
- `FE_SimplexP<2>(k)` with configurable polynomial degree `k` (`fe_degree >= 1`)

## Code Structure

- [CMakeLists.txt](/home/chenyuye/homework/pde_project/CMakeLists.txt): build configuration
- [src/main.cpp](/home/chenyuye/homework/pde_project/src/main.cpp): CLI parsing, solver setup, run entry
- [src/Wave2D.hpp](/home/chenyuye/homework/pde_project/src/Wave2D.hpp): solver class, enums, and interfaces
- [src/Wave2D.cpp](/home/chenyuye/homework/pde_project/src/Wave2D.cpp): matrix assembly, boundary handling, time stepping, output, error/energy evaluation
- `scripts/mesh-square.geo`: Gmsh geometry for unit square meshes
- `scripts/generate-meshes.sh`: batch mesh generation
- `scripts/plot_energy.py`: energy and relative drift plotting from CSV
- `scripts/plot_error.py`: L2/H1 error plotting from CSV
- `scripts/plot_3d_surface.py`: 3D surface visualization from VTU files
- `scripts/generate_plots.sh`: post-run plotting pipeline

## Build

Example (course/module environment):

```bash
cd /home/chenyuye/homework/pde_project
module load gcc-glibc dealii
mkdir -p build
cd build
cmake ..
make -j
```

Clean all run outputs:

```bash
cd /home/chenyuye/homework/pde_project/build
make clean-results
```

## Run

```bash
cd /home/chenyuye/homework/pde_project/build
./main [mesh.msh] [dt] [T] [output_every] [omega] [time_scheme] [mass_type] [fe_degree] [compute_error_each_step] [auto_plot]
```

Arguments:

1. `mesh.msh`: Gmsh mesh file path
2. `dt`: time step size
3. `T`: final time
4. `output_every`: output VTU every N steps (no output when it equal to 0)
5. `omega`: angular frequency (used for `eigenmode` exact solution)
6. `time_scheme`: `cd` (default) or `newmark`
7. `mass_type`: `lumped` (default) or `consistent`
8. `fe_degree`: polynomial degree (default `1`)
9. `compute_error_each_step`: `1` (default) to compute/write error every step, `0` to write only final error
10. `auto_plot`: `1` (default) to auto-run plotting scripts, `0` to skip auto plotting

Default-equivalent run:

```bash
./main ../mesh/mesh-square-h0.1.msh 0.01 2.0 1 4.442882938 cd lumped 1 1 1
```

## Output Files

All outputs are written under:

- `result/<run_config>/` (at the project root)

where `<run_config>` is automatically generated from input parameters:
mesh, time scheme, mass type, FE degree, `dt`, `T`, `output_every`, and `omega`.

Inside each run folder, outputs are tagged by method:

- `output-<mesh>-<method_tag>-<step>.vtu`
- `energy-<method_tag>.csv`
- `error-<method_tag>.csv`
- `energy-<method_tag>.png` (auto-generated after run)
- `error-<method_tag>.png` (auto-generated after run)

where `method_tag = <time_scheme>-<mass_type>-p<fe_degree>`.

`error-<method_tag>.csv` columns:
- `step,time,L2_error,H1_error`

## Notes

1. `cd` is explicit and requires a CFL-like step restriction (`dt/h` small enough).
2. `newmark` (`beta=0.25`, `gamma=0.5`) is unconditionally stable for linear problems but each step solves a linear system.
3. To discuss numerical dissipation and dispersion, compare `energy-*.csv` and `error-*.csv` across different `(time_scheme, mass_type, fe_degree, dt, h)` configurations.
4. After each run, `main` calls `scripts/generate_plots.sh <run_output_dir>`, and the shell script invokes `plot_energy.py`/`plot_error.py` for available CSV files.
5. You can also run plotting manually:
   `bash scripts/generate_plots.sh result/<run_config>`
6. For interactive display when running Python scripts manually, use `--show`, e.g. `python3 scripts/plot_energy.py path/to/energy.csv --show`.

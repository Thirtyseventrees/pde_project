# Project 2 — Solving the 2D Wave Equation
## Finite Element Solution Workflow

---

## 0. Problem Statement

We consider the **2D wave equation**

$$
\begin{cases}
\dfrac{\partial^2 u}{\partial t^2} - \Delta u = f & \text{in } \Omega, \\\\
u = g & \text{on } \partial\Omega, \\\\
u(x,0)=u_0 & \text{in } \Omega, \\\\
\dfrac{\partial u}{\partial t}(x,0)=u_1 & \text{in } \Omega.
\end{cases}
$$

Goal:

- Implement a **finite element solver**
- Select suitable **space and time discretization**
- Discuss numerical properties:
  - dissipation
  - dispersion
  - computational aspects

---

# 1. Treatment of Dirichlet Boundary Conditions

For non–homogeneous Dirichlet boundary conditions:

$$
u = g \quad \text{on } \partial\Omega
$$

Two common approaches exist.

---

## Method A — Direct Enforcement (Recommended)

At every time step:

- Boundary degrees of freedom (DOFs) are fixed as

$$
U_B(t)=g(t)
$$

- Modify the right-hand side accordingly.

✅ Simple  
✅ Widely used in implementations

---

## Method B — Homogenization

Introduce

$$
u = w + \tilde g
$$

where

$$
\tilde g = g \text{ on } \partial\Omega
$$

Then

$$
w|_{\partial\Omega}=0
$$

This simplifies theory but complicates implementation.

---

# 2. Weak Formulation

Multiply the PDE by a test function \(v\) and integrate:

$$
\int_\Omega u_{tt} v\,dx
+
\int_\Omega \nabla u \cdot \nabla v\,dx
=
\int_\Omega f v\,dx
$$

The boundary term vanishes because

$$
v=0 \quad \text{on } \partial\Omega
$$

---

# 3. Spatial Discretization (Finite Element Method)

---

## 3.1 Mesh and FE Space

- Triangular mesh: \( \mathcal{T}_h \)
- Choose FEM space:
  - P1 (linear)
  - P2 (quadratic)

Approximation:

$$
u_h(x,t)=\sum_{i=1}^{N} U_i(t)\phi_i(x)
$$

---

## 3.2 Matrix Formulation

After discretization:

$$
M\ddot U(t) + K U(t) = F(t)
$$

---

### Mass Matrix

$$
M_{ij}=\int_\Omega \phi_i \phi_j dx
$$

---

### Stiffness Matrix

$$
K_{ij}=\int_\Omega \nabla\phi_i\cdot\nabla\phi_j dx
$$

---

### Load Vector

$$
F_i(t)=\int_\Omega f(x,t)\phi_i dx
$$

---

## 3.3 Assembly Procedure

For each element \(T\):

1. Compute local matrices
   - \(M^T\)
   - \(K^T\)
2. Add contributions to global matrices

---

# 4. Initial Condition Projection

We need:

- displacement \(U(0)\)
- velocity \(\dot U(0)\)

---

## L2 Projection

Solve

$$
M U^0=b^0
$$

where

$$
b^0_i=\int_\Omega u_0\phi_i dx
$$

Similarly,

$$
M V^0=b^1
$$

---

## Nodal Interpolation (Simpler)

$$
U_i^0=u_0(x_i)
$$

$$
V_i^0=u_1(x_i)
$$

---

# 5. Time Discretization

The semi-discrete system:

$$
M\ddot U + KU = F
$$

---

# Method 1 — Explicit Central Difference

---

## Time Approximation

$$
\ddot U^n
\approx
\dfrac{U^{n+1}-2U^n+U^{n-1}}{\Delta t^2}
$$

Substitution gives:

$$
U^{n+1}
=
2U^n-U^{n-1}
+
\Delta t^2 M^{-1}(F^n-KU^n)
$$

---

## Mass Lumping

To avoid solving linear systems:

- Replace \(M\) with diagonal approximation
- Enables explicit update

---

## CFL Stability Condition

$$
\Delta t \le C \frac{h_{\min}}{c}
$$

---

## Numerical Properties

### Dispersion
- Phase velocity depends on wavelength
- Causes wave distortion

### Dissipation
- Very small for central difference
- Energy nearly conserved

---

# Method 2 — Newmark Scheme

Typical parameters:

$$
\gamma=\frac12,\quad
\beta=\frac14
$$

Solve each step:

$$
\left(
\frac{1}{\beta\Delta t^2}M + K
\right)
U^{n+1}
=
\text{RHS}
$$

### Properties

✅ Unconditionally stable  
✅ Controllable dissipation  
❌ Requires linear solver

---

# Method 3 — First Order Reformulation

Rewrite system:

$$
\dot U = V
$$

$$
M\dot V + KU = F
$$

Apply:

- Crank–Nicolson
- Backward Euler

---

# 6. Boundary Enforcement During Time Stepping

Partition DOFs:

$$
U=
\begin{bmatrix}
U_I \\\\
U_B
\end{bmatrix}
$$

Boundary values known:

$$
U_B=g(t)
$$

Implementation shortcut:

- Update solution
- Overwrite boundary DOFs

---

# 7. Verification and Error Analysis

---

## Manufactured Solution

Choose analytic solution:

$$
u_{exact}(x,y,t)
$$

Compute corresponding:

- forcing term
- boundary data
- initial conditions

---

## Error Measures

### L2 Error

$$
\|u_h-u\|_{L^2(\Omega)}
$$

### H1 Error

$$
\|\nabla(u_h-u)\|_{L^2(\Omega)}
$$

---

# 8. Numerical Dissipation and Dispersion Discussion

---

## Dispersion

Discrete waves satisfy modified dispersion relation:

- Different frequencies propagate at different speeds
- Causes phase error

Influenced by:
- mesh size
- time step
- element order

---

## Dissipation

Ideal wave equation conserves energy:

$$
E=
\frac12\dot U^T M\dot U
+
\frac12 U^T K U
$$

Some schemes introduce artificial damping.

---

# 9. Solver Algorithm Structure

---

## Step 1 — Mesh Generation

Create mesh and DOFs.

---

## Step 2 — Assemble Matrices

Compute:

- mass matrix \(M\)
- stiffness matrix \(K\)

---

## Step 3 — Initialization

Compute:

- \(U^0\)
- \(V^0\)

For explicit scheme:

$$
U^1
=
U^0
+\Delta t V^0
+\frac{\Delta t^2}{2}A^0
$$

---

## Step 4 — Time Loop

For each timestep:

1. Compute load vector
2. Update solution
3. Apply boundary conditions
4. Output results

---

## Step 5 — Postprocessing

- error computation
- energy monitoring
- visualization (VTK / ParaView)

---

# 10. Computational Complexity

Explicit scheme per step:

- sparse matrix-vector product
- vector operations

Cost:

$$
O(N)
$$

Suitable for large-scale simulations.

---
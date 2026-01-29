# Hydraulic Transient (Water Hammer) Simulation — Fortran (MOC)

This project implements a **Fortran-based numerical simulation of hydraulic transients (water hammer)** in a **pipeline system consisting of multiple pipes connected in series**. The model captures unsteady flow behavior caused by **time-dependent valve operations** at the downstream end, using a **Method of Characteristics (MOC)**–style formulation.

All calculations are performed using **SI units**.

---

## Features

- Supports up to **10 pipes connected in series**
- Each pipe may contain up to **100 computational nodes**
- **Wave speed adjustment** to eliminate interpolation errors
- **Darcy–Weisbach friction** included
- **Time-varying valve boundary condition**
- Valve opening defined by **relative opening (τ) vs time curve**
- **Parabolic interpolation** for intermediate valve positions
- Tracks **maximum and minimum piezometric head** at every node
- Suitable for **water hammer and transient pressure analysis**

> Array limits can be increased by modifying the `dimension` statements.

---

## Numerical Method Overview

- The governing equations for unsteady flow are solved using a **time-marching scheme**
- Interior nodes are updated using **compatibility equations along characteristic lines**
- Boundary conditions:
  - **Upstream**: constant-head reservoir
  - **Downstream**: valve with time-dependent opening
- Friction losses are incorporated directly into the characteristic equations

---

## Program Limits

| Parameter | Current Limit |
|---------|---------------|
| Number of pipes | 10 |
| Nodes per pipe | 100 |
| Valve curve points | 10 |

To increase limits, modify the following declarations:

```fortran
real , dimension(10,100) :: q,h,qp,hp,hmax,hmin
real , dimension(10)     :: ca,f,cf,ar,a,l,n,d,y,an,aunadj,ani

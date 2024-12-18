# Documentation

## Contents

* [How it works](#How-it-works)
* [Installation](installation.md)
    > [For C projects and CLI](installation.md#For-C-projects-and-CLI)
    > [For Python projects](installation.md#For-Python-projects)
* [Using dofulator](usage.md)
    > [CLI](usage.md#CLI)
    > [C](usage.md#C)
    > [Python](usage.md#Python)
    > [MDAnalysis](usage.md/MDAnalysis)


## How it works

For a set of atoms, $\mathcal{S}$, the instantaneous local kinetic temperature is given by
$$T_{\mathcal{S}} = \frac{\sum_{i\in\mathcal{S}} m_i\mathbf{v}_i\cdot\mathbf{v}_i}{k_B\sum_{i\in\mathcal{S}} d_i},$$
where $m_i$ is the atomic mass, $\mathbf{v}_i$ is the velocity, $k_B$ is Boltzmann's
constant, and $d_i$ is the degrees of freedom (DoF) associated with atom $i$.
For unconstrained 3D systems, $d_i$ is simply 3. However, when constraints are
introduced, and those constraints include atoms both inside and outside of
$\mathcal{S}$, it becomes non-trivial to calculate $\sum_{i\in\mathcal{S}}
d_i$.

It turns out that atomic DoF can be determined from the Jacobian which maps
modal velocities of a fragment (group of atoms connected by at least 1 bond, or
atoms which form a rigid body) to atomic velocities.
This library groups atoms into fragments based on provided connectivity
information, and then calculates the Jacobian for a given conformation and uses
it to determine $d_i$.
Optionally, $d_i$ can be broken down into separate $x$, $y$ and $z$ directional
contributions.
For details of the theory, see
[J. Chem. Theory Comput. 2024, 20, 23, 10615-10624](https://doi.org/10.1021/acs.jctc.4c00957).

For performance and compatibility reasons, dofulator is written as a C library,
but a python wrapper is also provided for convenience which includes an
[MDAnalysis](https://github.com/MDAnalysis/mdanalysis) plugin.
A basic CLI is also included, which is primarily intended for quickly calculating
atomic DoF of single molecules.

Note, currently the global center of mass constraint which results in a loss of
1 DoF in each direction over the entire system is currently not accounted for,
but the resultant error is negligible except in very small systems.

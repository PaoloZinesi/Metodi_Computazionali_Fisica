# Metodi Computazionali della Fisica ![ ](https://img.shields.io/badge/c++-00599d?style=for-the-badge&logo=cplusplus&logoColor=white)

Repository for the Bachelor course in Computational Methods of Physics (Metodi Computazionali della Fisica), year 2020/2021.
This repository contains all my laboratory assignments and the code presented as my final exam.

### Topics of the lab assignments
- **Lab1:** Grid definition and function discretization
- **Lab2:** Numerical derivatives and numerical integrals
- **Lab3:** Harmonic oscillator solution with Runge-Kutta method
- **Lab4:** Study of Fermi-Pasta-Ulam problem with Velocity Verlet method
- **Lab5:** 2D Poisson equation with iterative Jacobi method
- **Lab6:** 1D time-dependent Schrödinger equation with Crank-Nicolson method
- **Lab7:** Ising model with Metropolis algorithm

### Exam topic
The final exam is a solution of the 2D time-dependent Schrödinger equation describing a particle confined in a rectangular box.
This project employs the open-source library [Eigen](https://gitlab.com/libeigen/eigen/) to optimize matrix equations when matrices are sparse (which is typical in the numerical solution of PDEs).
- _v1_ contains the code written in a single main program
- _v2_ contains the code subdivided into multiple routines

Some sample images are shown below. Here, the evolving wave function is subject to a potential that emulates a double-slit. As expected, the results are consistent with Young's interference experiment.
<p align="center">
    <table border="0">
        <tr>
            <td> <img src="/Esame/images/young_high_1.jpg" width="480"> </td>
            <td> <img src="/Esame/images/young_high_2.jpg" width="480"> </td>
        </tr>
        <tr>
            <td> <img src="/Esame/images/young_high_3.jpg" width="480"> </td>
            <td> <img src="/Esame/images/young_high_4.jpg" width="480"> </td>
        </tr>
    </table>
</p>
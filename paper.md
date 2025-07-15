---
title: '<span style="font-family:qcr;">nudustc++</span>: C++ Code for Modeling Dust Nucleation and Destruction in Gaseous Systems'
tags:
  - astronomy
  - nucleation
  - C++
  - Supernovae
  - astrophysics
  - dust
  - dust grains
  - grain physics
  - physics
authors:
  - name: Sarah M. Stangl
    orcid: 0000-0001-5570-6666
    equal-contrib: true
    affiliation: "1, 2,3" 
  - name: Christopher M. Mauney
    orcid: 0000-0002-7827-2247
    equal-contrib: true 
    affiliation: "1,2,4"
affiliations:
 - name: Theoretical Design Division, Los Alamos National Laboratory, Los Alamos, NM, 87545, USA
   index: 1
 - name: Center for Theoretical Astrophysics, Los Alamos National Laboratory, Los Alamos, NM, 87545, USA
   index: 2
 - name: University of Oklahoma, Norman, OK, 73019, USA
   index: 3
 - name: High Performance Computing Division, Los Alamos National Laboratory, Los Alamos, NM, 87545, USA
   index: 4
date: 8 January 2024
bibliography: paper.bib
---

# Summary

We introduce <span style="font-family:qcr;">nudustc++</span>, a **nu**cleating **dust** code in **C++** modeling dust grain formation, growth, and erosion in gaseous systems. <span style="font-family:qcr;">nudustc++</span> is a highly parallelizable set of code and tools for solving a system of nonlinear ordinary differential equations describing dust nucleation, growth, and erosion for user-specified grain species. It leverages OpenMP and MPI to optimize threading and distribution on available CPUs.

# Statement of need

Understanding interstellar dust is crucial for astronomical observations [@Draine2003], offering key insights into stellar processes. These grains absorb electromagnetic radiation, re-emitting it at longer wavelengths, leading to extinction and a spectral shift towards redder wavelengths. The size and composition of dust introduces variability in opacities and distortion of incident light, resulting in molecular lines and altering the resulting data. Dust formed in asymptotic giant branch (AGB) stars, on pre-existing grains in the interstellar medium (ISM), and within the expanding, cooling ejecta of core collapse supernovae explosions (CCSNe); these grains preserve important information about the nucleosynthetic processes within their host environment, locking up their unique isotopic signatures. However, despite their importance, the quantity, composition, and size distribution of dust formed in supernovae and deposited into the ISM remain poorly constrained. 

Models of formation, growth, and weathering of dust are necessary to understand the origin and characteristics of dust in the interstellar medium, shedding light on where and what dust are possible in these environments. Specifically, modeling the formation and survivability of dust in CCSNe produce a population of dust grains that can be compared with observations, allowing verification of our current understanding of physical models: ISM dust origin, chemical networks, nucleation models, hydrodynamics, supernovae engines, progenitor structure, erosion physics, and stellar compositions. 

This project originated from the need to track dust nucleation and destruction in Core-Collapse Supernovae Explosions. It addresses the lack of sub-grid physics associated with phase transitions, where hydrodynamical code's timesteps are an order of magnitude larger than those needed to capture gas vapor physics. A smaller more refined grid with chemical networks and smaller timesteps is needed. The code is structured to intake any hydro-dynamical temperature-density trajectory with vapor compositions. This allows <span style="font-family:qcr;">nudustc++</span> to track dust in a large range of environments: planetary atmospheres, nebulae, hydro-aerosoal formation, explosions, etc. Additionally, if a time series hydro-dynamical profile is unavailable but a dust size distribution and a profile snapshot is, <span style="font-family:qcr;">nudustc++</span> can calculate the evolution and survivability of the dust. The applications of <span style="font-family:qcr;">nudustc++</span> extend beyond Supernovae and Astronomy to include any model with thermodynamic and statistical physics. 

# State Of The Field

In order to gain a deeper understanding of the origin and characteristics of dust in the interstellar medium, it is imperative to develop models that include the nucleation, chemistry, growth, and erosion of dust. Current methods of calculating dust formation and survival include Classical Nucleation Theory (CNT) and Kinetic Nucleation Theory (KNT).

CNT treats grain formation as a barrier crossing problem. As atoms stick to a cluster, the free energy increases. After reaching a critical size, the free energy decreases as atoms are added. It tracks the nucleation rate by assuming a steady state between monomer attachment and detachment. However, it neglects chemical reactions of formation, destruction, growth by coagulation, and treats the grains as bulk materials. Due to these simplified assumptions, CNT is widely used, but is increasingly less so due to these limitations.  @kozasa1987 and @Bianchi2007 used CNT to model dust grain formation in SN 1987A and SN 1987A-like Supernovae. @Todini2001 used CNT to calculate dust formation in Core Collapse Supernovae Explosions. More recently @paquette2023 used CNT to model dust formation in the outflows of AGB Stars. 

KNT tracks the number densities of clusters with more than two atoms, treating them as particles that grow as spheres. It uses size dependent grain properties in place of bulk material properties. Grains grow by accreting atoms (condensation) through kinetic theory. Grains lose atoms through erosion (destruction) using the principle of detailed balance. While KNT does not assume a steady state between condensation and destruction, it still doesn't take into account the chemical reactions undergone or growth through coagulation. @Nozawa_2003 used KNT to model the nucleation and growth of dust in early Population III star supernovae. @Nozawa2006 included destruction to model the effects of high velocity shock waves on dust. @fallest2012 predicts dust mass yields in CCSNe using KNT.

Other publicly available dust codes include <span style="font-family:qcr;">sndust</span>, <span style="font-family:qcr;">starchem</span>, <span style="font-family:qcr;">DustPy</span> and <span style="font-family:qcr;">astrochem</span>. <span style="font-family:qcr;">sndust</span> is described in @brooker2021 and is a post processing dust nucleation code using KNT. <span style="font-family:qcr;">starchem</span> tracks the chemical network in stellar enviornments [@starchem]. <span style="font-family:qcr;">DustPy</span> is a python package used to model dust evolution in protoplanetary disks [@dustpy]. <span style="font-family:qcr;">astrochem</span> computes the chemical abundances and includes gas-dust interactions in astronomical environments such as the interstellar medium, diffuse clouds, and protostars [@astrochem]. 

# Design Principles and Salient Features

<span style="font-family:qcr;">nudustc++</span> is designed as a flexible, multi-use post-processing code. It is designed to take user supplied data in ascii, binary, or text format, including initial conditions, environment variables, and dust formation networks, and calculate a solution vector. We use an initial composition and interpolated user input data to construct a rate-of-change vector supplied to an implicit integrator in each cell. The user can easily change and modify the interpolator and integrator. Three main configuration paths are currently implemented: destruction only, destruction with nucleation and growth, and nucleation and growth. This allows the reduction of total computations by removing from the solution vector parameters modeling quantities not needed for the calculation path. Output data can be written to ascci, binary, or text files.

Because of the post processing nature of <span style="font-family:qcr;">nudustc++</span>, modeling grain nucleation and destruction is possible in a large range of physical environments. The user provides the hydrodynamical trajectory file and vapor compositions which can describe Supernovae explosions to planetary atmospheres to interstellar gas clouds and any temperature-density profile extended beyond astrophysics. 

Because <span style="font-family:qcr;">nudustc++</span> is a post processing code where data is read in separately for each cell with no data flow between cells, there is no shared memory and each cell can be computed independently. This results in an embarrassingly parallel process, allowing for simultaneous computation of each cell and thereby reducing runtime. 

# Performance

The Figure below shows a comparison between the python <span style="font-family:qcr;">sndust</span> using an implicit integrator versus <span style="font-family:qcr;">nudustc++</span> using an explicit integrator. At shorter solution lengths, the implicit integrator has a lot of overhead leading to an increased run time despite the shorter length. The explicit integrator in <span style="font-family:qcr;">C++</span> out performs the implicit integrator even at longer solution length. Overall, <span style="font-family:qcr;">nudustc++</span> outperforms the python version. It is also embarrassingly parallel. Utilizing OpenMP and MPI, the user can run many cells in parallel. 

# Figures

![Performance Comparison of <span style="font-family:qcr;">nudustc++</span> and implicit and explicit integrators.\label{fig:performance}](int_logtime_cycles.png)


# Acknowledgements

This work was supported by the U.S. Department of Energy through the Los Alamos National Laboratory (LANL). LANL is operated by Triad National Security, LLC, for the National Nuclear Security Administration of U.S. Department of Energy (Contract No. 89233218CNA000001). This research used resources provided by the Darwin testbed at LANL which is funded by the Computational Systems and Software Environments subprogram of LANL’s Advanced Simulation and Computing program (NNSA/DOE). Work supported by the Laboratory Directed Research and Development program of Los Alamos National Laboratory, under project number  20240039DR. This work is approved for unlimited release with report number LA-UR-24-20630.

# References

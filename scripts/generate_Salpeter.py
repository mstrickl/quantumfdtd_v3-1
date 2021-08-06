#!/usr/bin/env python3

import os
import numpy as np

#####################################################################################################

def salp(fname, k, v, b, N, A):

    f = open(fname, 'w')

    for ikx in range(1,2*N+2):
        for iky in range(1,2*N+2):
            for ikz in range(1,2*N+2):

                kx = ikx-N
                ky = iky-N
                kz = ikz-N

                dx = kx+.5
                dy = ky+.5
                dz = kz+.5

                r = A*np.sqrt(dx*dx+dy*dy+dz*dz)

                salp = -k/r - v*np.exp(-b*r)/r

                f.write(f'{kx} {ky} {kz} {salp} 0.0\n')

    f.close()

#####################################################################################################

def gen_config(k, v, b, N=64, A=0.301, MASS=0.3):

    name=f'N{N}_A{A}_k{k}_v{v}_b{b}'
    os.makedirs(f'{name}/data/snapshot')
    os.makedirs(f'{name}/input')
    salp(f'{name}/input/pot.txt', k, v, b, N, A)
    f_param=f'{name}/input/params.txt'

    param_file = f"""// standard deviation of initial wavefunc noise 
SIG		1.0

// number of spatial grid points; should be divisible by the number of computational nodes
NUM 		{N}

// spatial grid spacing // default unit is 1/GeV
A               {A}	

// temporal grid spacing - should be <= A*A/3
EPS		.002

// convergence tolerance for ground state energy
// negative value for disabling the tolerance condition,
// so that the STEPS parameter will be the actual # of steps to take
TOLERANCE	-1.

// maximum # of steps to take
// If TOLERANCE is negative, this will be the actual # of steps to take
//STEPS   	9000
STEPS   	40000

// how many steps before recording measurables
UPDATE		10

// how many steps before dumping all variables to file; "taking a snapshot"
SNAPUPDATE	50

// set to one to dump debugging files containing snapshot "slices"
DUMPSNAPS	4

// Number of steps before recording a snap
SNAPDUMP        500

// set to one to dump the full 3d wavefncs to disk at the end of the run
SAVEWAVEFNCS	1

// set to one to dump the decay table time -> energy
SAVEDECAY       1

// set to one to use the relativistic version of the kintetic term: 20/03/2018
KINTERM         3

//
// potential to simulate
//
//    0  No potential (V=0); solutions to infinite-depth 3d well due to BCs
//    1  Cubic well
//    2  Radial Coulomb
//    3  Elliptical Coulomb
//    4  3d harmonic oscillator
//    5  Complex 3d hoarmonic oscillator
//    6  Cornell
//    90 Read potential.txt file, with rows (x,y,z,ReV,ImV)
//    91 Read potential.txt file, with rows (r2,ReV,ImV)

POTENTIAL    90

// External potential file (POTENTIAL fro 90 to 99 and 190 to 199)
EXTPOT       input/pot.txt

// set to one to save potential to file
SAVEPOT      1

// if set, and POTENTIAL=90, the potential is a chi^2 adjustment to V(r)=A+B/r+s*r for
// (r/A)>POTCRITR
POTCRITR  10000

// if set, and POTENTIAL=90, the potential is flat for (r/A)>POTFLATR.
POTFLATR  10000	

//
// initial condition to use
//
//    0         read initial condition from wavefunction_0_#.dat files in the data directory,
//              necessary files can be generated at end of a run by turning SAVEWAVEFNCS on
//    1 	random gaussian with std dev SIG
//    2		coulomb-like
//    3		constant of 0.1 in interior 
//    4		boolean test grid; mod(i%2)mod(j%2)mod(k%2)
//
INITCONDTYPE	2

//
// initial symmetry constraint
//
//	0	None
//	1	Symmetric about z-axis
//	2	Antisymmetric about z-axis
//	3	Symmetric about y-axis
//	4	Antisymmetric about y-axis
//
INITSYMMETRY	0

//
// Physical parameters used in potentials
//

// Reduced mass; charmonium reduced mass is 1.29/2; bottomonium reduced mass is 4.7/2
MASS		{MASS}

// string tension in units GeV^2
SIGMA		0.223
"""
    with open(f_param,'w') as f: f.write(param_file)


#####################################################################################################

gen_config(0.5, 0.5, 0.5)

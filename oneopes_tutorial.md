# Complete OneOPES Tutorial for Host-Guest Systems

## Table of Contents
1. [System Preparation](#1-system-preparation)
2. [Initial PLUMED Setup for Unbiased Run](#2-initial-plumed-setup-for-unbiased-run)
3. [Running Unbiased Simulation](#3-running-unbiased-simulation)
4. [Setting Up OneOPES Replicas](#4-setting-up-oneopes-replicas)
5. [Running OneOPES Simulation](#5-running-oneopes-simulation)

## 1. System Preparation

1. Calculate atomic Charges
2. Parameterize System
3. Obtain initial binding pose
4. Generate topology files
5. Add water model and ions
6. Run unbiased simulations

    1.6.1. Energy minimization

    1.6.2. NVT equilibration

    1.6.3. NPT equilibration

    1.6.4. Generate .tpr file for production

## 2. Initial PLUMED Setup for Unbiased Run

### 2.1 Define Groups
First, identify correct atom numbers from your .gro file:

```plumed
# Host, ligand and water groups
HOST: GROUP ATOMS=29-122         # Host atoms
LIGC: GROUP ATOMS=1-11           # Ligand heavy atoms
WO: GROUP ATOMS=122-4235:3       # Water oxygen atoms

# Selected ligand atoms for orientation and water coordination
l1: GROUP ATOMS=1                # Ligand atom 1
l2: GROUP ATOMS=8                # Ligand atom 2
l3: GROUP ATOMS=7                # Ligand atom 3
l4: GROUP ATOMS=10               # Ligand atom 4

# System treatment
WHOLEMOLECULES ENTITY0=HOST
FIT_TO_TEMPLATE STRIDE=1 REFERENCE=host_template.pdb TYPE=OPTIMAL
lig: CENTER ATOMS=LIGC
```

### 2.2 Set Up Virtual Atoms and Basic CVs
Place virtual atoms along the Z-axis of the host
```plumed
#HIDDEN
# Host, ligand and water groups
HOST: GROUP ATOMS=28-171         # Host atoms
LIGC: GROUP ATOMS=1-11           # Ligand heavy atoms
WO: GROUP ATOMS=172-4234:3       # Water oxygen atoms

# Selected ligand atoms for orientation and water coordination
l1: GROUP ATOMS=1                # Ligand atom 1
l2: GROUP ATOMS=8                # Ligand atom 2
l3: GROUP ATOMS=7                # Ligand atom 3
l4: GROUP ATOMS=10               # Ligand atom 4

# System treatment
WHOLEMOLECULES ENTITY0=HOST
FIT_TO_TEMPLATE STRIDE=1 REFERENCE=host_template.pdb TYPE=OPTIMAL
lig: CENTER ATOMS=LIGC
#ENDHIDDEN
# Virtual atoms defining funnel axis
v1: FIXEDATOM AT=0.000,0.000,0.000
v2: FIXEDATOM AT=0.000,0.000,0.250
v3: FIXEDATOM AT=0.000,0.000,0.500
v4: FIXEDATOM AT=0.000,0.000,0.750
v5: FIXEDATOM AT=0.000,0.000,1.000
v6: FIXEDATOM AT=0.000,0.000,-0.250
v7: FIXEDATOM AT=0.000,0.000,-0.500
v8: FIXEDATOM AT=0.000,0.000,-0.750
v9: FIXEDATOM AT=0.000,0.000,-1.000

# Basic position variables
cyl: DISTANCE ATOMS=v1,lig COMPONENTS
radius: MATHEVAL ARG=cyl.x,cyl.y FUNC=sqrt(x*x+y*y) PERIODIC=NO
```

### 2.3 Define Water Coordination
Define water coordination between water oxygen atoms and the selected ligand atoms and virtual atoms. 
```plumed
#HIDDEN
# Host, ligand and water groups
HOST: GROUP ATOMS=28-171         # Host atoms
LIGC: GROUP ATOMS=1-11           # Ligand heavy atoms
WO: GROUP ATOMS=172-4234:3       # Water oxygen atoms

# Selected ligand atoms for orientation and water coordination
l1: GROUP ATOMS=1                # Ligand atom 1
l2: GROUP ATOMS=8                # Ligand atom 2
l3: GROUP ATOMS=7                # Ligand atom 3
l4: GROUP ATOMS=10               # Ligand atom 4

# System treatment
WHOLEMOLECULES ENTITY0=HOST
FIT_TO_TEMPLATE STRIDE=1 REFERENCE=host_template.pdb TYPE=OPTIMAL
lig: CENTER ATOMS=LIGC

# Virtual atoms defining funnel axis
v1: FIXEDATOM AT=0.000,0.000,0.000
v2: FIXEDATOM AT=0.000,0.000,0.250
v3: FIXEDATOM AT=0.000,0.000,0.500
v4: FIXEDATOM AT=0.000,0.000,0.750
v5: FIXEDATOM AT=0.000,0.000,1.000
v6: FIXEDATOM AT=0.000,0.000,-0.250
v7: FIXEDATOM AT=0.000,0.000,-0.500
v8: FIXEDATOM AT=0.000,0.000,-0.750
v9: FIXEDATOM AT=0.000,0.000,-1.000

# Basic position variables
cyl: DISTANCE ATOMS=v1,lig COMPONENTS
radius: MATHEVAL ARG=cyl.x,cyl.y FUNC=sqrt(x*x+y*y) PERIODIC=NO
#ENDHIDDEN
# Ligand-water coordination
L1: COORDINATION GROUPA=l1 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
L2: COORDINATION GROUPA=l2 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
L3: COORDINATION GROUPA=l3 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
L4: COORDINATION GROUPA=l4 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20

# Virtual atom-water coordination
V1: COORDINATION GROUPA=v1 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V2: COORDINATION GROUPA=v2 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V3: COORDINATION GROUPA=v3 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V4: COORDINATION GROUPA=v4 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V5: COORDINATION GROUPA=v5 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V6: COORDINATION GROUPA=v6 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V7: COORDINATION GROUPA=v7 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V8: COORDINATION GROUPA=v8 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V9: COORDINATION GROUPA=v9 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
```

### 2.4 Define Funnel, Walls, and Additional CVs
```plumed
#HIDDEN
# Host, ligand and water groups
HOST: GROUP ATOMS=28-171         # Host atoms
LIGC: GROUP ATOMS=1-11           # Ligand heavy atoms
WO: GROUP ATOMS=172-4234:3       # Water oxygen atoms

# Selected ligand atoms for orientation and water coordination
l1: GROUP ATOMS=1                # Ligand atom 1
l2: GROUP ATOMS=8                # Ligand atom 2
l3: GROUP ATOMS=7                # Ligand atom 3
l4: GROUP ATOMS=10               # Ligand atom 4

# System treatment
WHOLEMOLECULES ENTITY0=HOST
FIT_TO_TEMPLATE STRIDE=1 REFERENCE=host_template.pdb TYPE=OPTIMAL
lig: CENTER ATOMS=LIGC

# Virtual atoms defining funnel axis
v1: FIXEDATOM AT=0.000,0.000,0.000
v2: FIXEDATOM AT=0.000,0.000,0.250
v3: FIXEDATOM AT=0.000,0.000,0.500
v4: FIXEDATOM AT=0.000,0.000,0.750
v5: FIXEDATOM AT=0.000,0.000,1.000
v6: FIXEDATOM AT=0.000,0.000,-0.250
v7: FIXEDATOM AT=0.000,0.000,-0.500
v8: FIXEDATOM AT=0.000,0.000,-0.750
v9: FIXEDATOM AT=0.000,0.000,-1.000

# Basic position variables
cyl: DISTANCE ATOMS=v1,lig COMPONENTS
radius: MATHEVAL ARG=cyl.x,cyl.y FUNC=sqrt(x*x+y*y) PERIODIC=NO

# Ligand-water coordination
L1: COORDINATION GROUPA=l1 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
L2: COORDINATION GROUPA=l2 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
L3: COORDINATION GROUPA=l3 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
L4: COORDINATION GROUPA=l4 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20

# Virtual atom-water coordination
V1: COORDINATION GROUPA=v1 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V2: COORDINATION GROUPA=v2 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V3: COORDINATION GROUPA=v3 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V4: COORDINATION GROUPA=v4 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V5: COORDINATION GROUPA=v5 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V6: COORDINATION GROUPA=v6 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V7: COORDINATION GROUPA=v7 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V8: COORDINATION GROUPA=v8 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V9: COORDINATION GROUPA=v9 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
#ENDHIDDEN
# Double funnel for symmetric hosts (e.g., CB8)
funnel: MATHEVAL ARG=radius,cyl.z VAR=r,z FUNC=((r+1.5*(-2.0+z))*step(-z+0.5)+(r-0.2)*step(z-0.5))*step(z-0.0)+((r-1.5*(2.0+z))*step(z+0.5)+(r-0.2)*step(-z-0.5))*step(-z+0.0) PERIODIC=NO

# Walls
UPPER_WALLS AT=0 ARG=funnel KAPPA=2000.0 LABEL=funnelwall
UPPER_WALLS AT=1.5 ARG=cyl.z KAPPA=2000.0 LABEL=upper_wall
LOWER_WALLS AT=-1.5 ARG=cyl.z KAPPA=2000.0 LABEL=lower_wall

# Angle and energy
ang: ANGLE ATOMS=v3,v5,1,10
cosang: MATHEVAL ARG=ang FUNC=cos(x) PERIODIC=NO
ene: ENERGY

# Print all variables for sigma calculation
PRINT ARG=cyl.z,radius,cosang,L1,L2,L3,L4,V1,V2,V3,V4,V5,V6,V7 STRIDE=50 FILE=COLVAR FMT=%8.4f
```

## 3. Running Unbiased Simulation

1. Run 1 ns simulation:
```bash
gmx_mpi mdrun -deffnm topol -plumed plumed.dat -nsteps 500000
```

2. Calculate standard deviations of CVs from COLVAR file
3. These values will be used as SIGMA values in the biased simulations

## 4. Setting Up OneOPES Replicas

Create 8 replica directories and prepare a plumed.dat file for each:

### Replica 0 (Base)
```plumed
# Include all group definitions and CVs from unbiased setup
#HIDDEN
# Host, ligand and water groups
HOST: GROUP ATOMS=28-171         # Host atoms
LIGC: GROUP ATOMS=1-11           # Ligand heavy atoms
WO: GROUP ATOMS=172-4234:3       # Water oxygen atoms

# Selected ligand atoms for orientation and water coordination
l1: GROUP ATOMS=1                # Ligand atom 1
l2: GROUP ATOMS=8                # Ligand atom 2
l3: GROUP ATOMS=7                # Ligand atom 3
l4: GROUP ATOMS=10               # Ligand atom 4

# System treatment
WHOLEMOLECULES ENTITY0=HOST
FIT_TO_TEMPLATE STRIDE=1 REFERENCE=host_template.pdb TYPE=OPTIMAL
lig: CENTER ATOMS=LIGC

# Virtual atoms defining funnel axis
v1: FIXEDATOM AT=0.000,0.000,0.000
v2: FIXEDATOM AT=0.000,0.000,0.250
v3: FIXEDATOM AT=0.000,0.000,0.500
v4: FIXEDATOM AT=0.000,0.000,0.750
v5: FIXEDATOM AT=0.000,0.000,1.000
v6: FIXEDATOM AT=0.000,0.000,-0.250
v7: FIXEDATOM AT=0.000,0.000,-0.500
v8: FIXEDATOM AT=0.000,0.000,-0.750
v9: FIXEDATOM AT=0.000,0.000,-1.000

# Basic position variables
cyl: DISTANCE ATOMS=v1,lig COMPONENTS
radius: MATHEVAL ARG=cyl.x,cyl.y FUNC=sqrt(x*x+y*y) PERIODIC=NO

# Ligand-water coordination
L1: COORDINATION GROUPA=l1 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
L2: COORDINATION GROUPA=l2 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
L3: COORDINATION GROUPA=l3 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
L4: COORDINATION GROUPA=l4 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20

# Virtual atom-water coordination
V1: COORDINATION GROUPA=v1 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V2: COORDINATION GROUPA=v2 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V3: COORDINATION GROUPA=v3 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V4: COORDINATION GROUPA=v4 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V5: COORDINATION GROUPA=v5 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V6: COORDINATION GROUPA=v6 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V7: COORDINATION GROUPA=v7 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V8: COORDINATION GROUPA=v8 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V9: COORDINATION GROUPA=v9 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
# Double funnel for symmetric hosts (e.g., CB8)
funnel: MATHEVAL ARG=radius,cyl.z VAR=r,z FUNC=((r+1.5*(-2.0+z))*step(-z+0.5)+(r-0.2)*step(z-0.5))*step(z-0.0)+((r-1.5*(2.0+z))*step(z+0.5)+(r-0.2)*step(-z-0.5))*step(-z+0.0) PERIODIC=NO

# Walls
UPPER_WALLS AT=0 ARG=funnel KAPPA=2000.0 LABEL=funnelwall
UPPER_WALLS AT=1.5 ARG=cyl.z KAPPA=2000.0 LABEL=upper_wall
LOWER_WALLS AT=-1.5 ARG=cyl.z KAPPA=2000.0 LABEL=lower_wall

# Angle and energy
ang: ANGLE ATOMS=v3,v5,1,10
cosang: MATHEVAL ARG=ang FUNC=cos(x) PERIODIC=NO
ene: ENERGY

#ENDHIDDEN
# Then add the main CVs:
# --- (4) OPES  ---
OPES_METAD_EXPLORE ...
    LABEL=opes
    ARG=cyl.z,cosang
    SIGMA=0.1,0.2    # Use values from unbiased run
    FILE=Kernels.data
    STATE_RFILE=compressed.Kernels.data
    STATE_WFILE=compressed.Kernels.data
    PACE=10000
    BARRIER=100
... OPES_METAD_EXPLORE

# Then update the print statement
PRINT ARG=opes.bias,cyl.z,radius,funnelwall.bias,upper_wall.bias,lower_wall.bias,cosang,L1,L2,L3,L4,V1,V2,V3,V4,V5,V6,V7,V8,V9 STRIDE=50 FILE=COLVAR FMT=%8.4f
```

### Replica 1
Add CV to previous replica:
```plumed
#HIDDEN
# Host, ligand and water groups
HOST: GROUP ATOMS=28-171         # Host atoms
LIGC: GROUP ATOMS=1-11           # Ligand heavy atoms
WO: GROUP ATOMS=172-4234:3       # Water oxygen atoms

# Selected ligand atoms for orientation and water coordination
l1: GROUP ATOMS=1                # Ligand atom 1
l2: GROUP ATOMS=8                # Ligand atom 2
l3: GROUP ATOMS=7                # Ligand atom 3
l4: GROUP ATOMS=10               # Ligand atom 4

# System treatment
WHOLEMOLECULES ENTITY0=HOST
FIT_TO_TEMPLATE STRIDE=1 REFERENCE=host_template.pdb TYPE=OPTIMAL
lig: CENTER ATOMS=LIGC

# Virtual atoms defining funnel axis
v1: FIXEDATOM AT=0.000,0.000,0.000
v2: FIXEDATOM AT=0.000,0.000,0.250
v3: FIXEDATOM AT=0.000,0.000,0.500
v4: FIXEDATOM AT=0.000,0.000,0.750
v5: FIXEDATOM AT=0.000,0.000,1.000
v6: FIXEDATOM AT=0.000,0.000,-0.250
v7: FIXEDATOM AT=0.000,0.000,-0.500
v8: FIXEDATOM AT=0.000,0.000,-0.750
v9: FIXEDATOM AT=0.000,0.000,-1.000

# Basic position variables
cyl: DISTANCE ATOMS=v1,lig COMPONENTS
radius: MATHEVAL ARG=cyl.x,cyl.y FUNC=sqrt(x*x+y*y) PERIODIC=NO

# Ligand-water coordination
L1: COORDINATION GROUPA=l1 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
L2: COORDINATION GROUPA=l2 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
L3: COORDINATION GROUPA=l3 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
L4: COORDINATION GROUPA=l4 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20

# Virtual atom-water coordination
V1: COORDINATION GROUPA=v1 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V2: COORDINATION GROUPA=v2 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V3: COORDINATION GROUPA=v3 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V4: COORDINATION GROUPA=v4 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V5: COORDINATION GROUPA=v5 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V6: COORDINATION GROUPA=v6 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V7: COORDINATION GROUPA=v7 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V8: COORDINATION GROUPA=v8 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V9: COORDINATION GROUPA=v9 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
# Double funnel for symmetric hosts (e.g., CB8)
funnel: MATHEVAL ARG=radius,cyl.z VAR=r,z FUNC=((r+1.5*(-2.0+z))*step(-z+0.5)+(r-0.2)*step(z-0.5))*step(z-0.0)+((r-1.5*(2.0+z))*step(z+0.5)+(r-0.2)*step(-z-0.5))*step(-z+0.0) PERIODIC=NO

# Walls
UPPER_WALLS AT=0 ARG=funnel KAPPA=2000.0 LABEL=funnelwall
UPPER_WALLS AT=1.5 ARG=cyl.z KAPPA=2000.0 LABEL=upper_wall
LOWER_WALLS AT=-1.5 ARG=cyl.z KAPPA=2000.0 LABEL=lower_wall

# Angle and energy
ang: ANGLE ATOMS=v3,v5,1,10
cosang: MATHEVAL ARG=ang FUNC=cos(x) PERIODIC=NO
ene: ENERGY

# --- (4) OPES  ---
OPES_METAD_EXPLORE ...
    LABEL=opes
    ARG=cyl.z,cosang
    SIGMA=0.1,0.2    # Use values from unbiased run
    FILE=Kernels.data
    STATE_RFILE=compressed.Kernels.data
    STATE_WFILE=compressed.Kernels.data
    PACE=10000
    BARRIER=100
... OPES_METAD_EXPLORE
#ENDHIDDEN
OPES_METAD_EXPLORE ...
    LABEL=opese1
    ARG=L4
    SIGMA=0.3    # From unbiased run
    FILE=Kernels1.data
    STATE_RFILE=compressed_Kernels1.data
    STATE_WFILE=compressed.Kernels1.data
    PACE=20000
    BARRIER=3
... OPES_METAD_EXPLORE

# Update print statement
PRINT ARG=opes.bias,cyl.z,radius,funnelwall.bias,upper_wall.bias,lower_wall.bias,cosang,L1,L2,L3,L4,V1,V2,V3,V4,V5,V6,V7,V8,V9,opese1.bias STRIDE=50 FILE=COLVAR FMT=%8.4f
```
### Replica 2
Add CV to previous replica:
```plumed
#HIDDEN
# Host, ligand and water groups
HOST: GROUP ATOMS=28-171         # Host atoms
LIGC: GROUP ATOMS=1-11           # Ligand heavy atoms
WO: GROUP ATOMS=172-4234:3       # Water oxygen atoms

# Selected ligand atoms for orientation and water coordination
l1: GROUP ATOMS=1                # Ligand atom 1
l2: GROUP ATOMS=8                # Ligand atom 2
l3: GROUP ATOMS=7                # Ligand atom 3
l4: GROUP ATOMS=10               # Ligand atom 4

# System treatment
WHOLEMOLECULES ENTITY0=HOST
FIT_TO_TEMPLATE STRIDE=1 REFERENCE=host_template.pdb TYPE=OPTIMAL
lig: CENTER ATOMS=LIGC

# Virtual atoms defining funnel axis
v1: FIXEDATOM AT=0.000,0.000,0.000
v2: FIXEDATOM AT=0.000,0.000,0.250
v3: FIXEDATOM AT=0.000,0.000,0.500
v4: FIXEDATOM AT=0.000,0.000,0.750
v5: FIXEDATOM AT=0.000,0.000,1.000
v6: FIXEDATOM AT=0.000,0.000,-0.250
v7: FIXEDATOM AT=0.000,0.000,-0.500
v8: FIXEDATOM AT=0.000,0.000,-0.750
v9: FIXEDATOM AT=0.000,0.000,-1.000

# Basic position variables
cyl: DISTANCE ATOMS=v1,lig COMPONENTS
radius: MATHEVAL ARG=cyl.x,cyl.y FUNC=sqrt(x*x+y*y) PERIODIC=NO

# Ligand-water coordination
L1: COORDINATION GROUPA=l1 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
L2: COORDINATION GROUPA=l2 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
L3: COORDINATION GROUPA=l3 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
L4: COORDINATION GROUPA=l4 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20

# Virtual atom-water coordination
V1: COORDINATION GROUPA=v1 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V2: COORDINATION GROUPA=v2 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V3: COORDINATION GROUPA=v3 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V4: COORDINATION GROUPA=v4 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V5: COORDINATION GROUPA=v5 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V6: COORDINATION GROUPA=v6 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V7: COORDINATION GROUPA=v7 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V8: COORDINATION GROUPA=v8 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V9: COORDINATION GROUPA=v9 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
# Double funnel for symmetric hosts (e.g., CB8)
funnel: MATHEVAL ARG=radius,cyl.z VAR=r,z FUNC=((r+1.5*(-2.0+z))*step(-z+0.5)+(r-0.2)*step(z-0.5))*step(z-0.0)+((r-1.5*(2.0+z))*step(z+0.5)+(r-0.2)*step(-z-0.5))*step(-z+0.0) PERIODIC=NO

# Walls
UPPER_WALLS AT=0 ARG=funnel KAPPA=2000.0 LABEL=funnelwall
UPPER_WALLS AT=1.5 ARG=cyl.z KAPPA=2000.0 LABEL=upper_wall
LOWER_WALLS AT=-1.5 ARG=cyl.z KAPPA=2000.0 LABEL=lower_wall

# Angle and energy
ang: ANGLE ATOMS=v3,v5,1,10
cosang: MATHEVAL ARG=ang FUNC=cos(x) PERIODIC=NO
ene: ENERGY

# --- (4) OPES  ---
OPES_METAD_EXPLORE ...
    LABEL=opes
    ARG=cyl.z,cosang
    SIGMA=0.1,0.2    # Use values from unbiased run
    FILE=Kernels.data
    STATE_RFILE=compressed.Kernels.data
    STATE_WFILE=compressed.Kernels.data
    PACE=10000
    BARRIER=100
... OPES_METAD_EXPLORE
OPES_METAD_EXPLORE ...
    LABEL=opese1
    ARG=L4
    SIGMA=0.3    # From unbiased run
    FILE=Kernels1.data
    STATE_RFILE=compressed_Kernels1.data
    STATE_WFILE=compressed.Kernels1.data
    PACE=20000
    BARRIER=3
... OPES_METAD_EXPLORE
#ENDHIDDEN
OPES_METAD_EXPLORE ...
    LABEL=opese2
    ARG=V1
    SIGMA=0.3   # From unbiased run
    FILE=Kernels2.data
    STATE_RFILE=compressed_Kernels2.data
    STATE_WFILE=compressed.Kernels2.data
    PACE=20000
    BARRIER=3
... OPES_METAD_EXPLORE

# Update print statement
PRINT ARG=opes.bias,cyl.z,radius,funnelwall.bias,upper_wall.bias,lower_wall.bias,cosang,L1,L2,L3,L4,V1,V2,V3,V4,V5,V6,V7,V8,V9,opese1.bias,opese2.bias STRIDE=50 FILE=COLVAR FMT=%8.4f
```

### Replica 3
Add CV to previous replica:
```plumed
#HIDDEN
# Host, ligand and water groups
HOST: GROUP ATOMS=28-171         # Host atoms
LIGC: GROUP ATOMS=1-11           # Ligand heavy atoms
WO: GROUP ATOMS=172-4234:3       # Water oxygen atoms

# Selected ligand atoms for orientation and water coordination
l1: GROUP ATOMS=1                # Ligand atom 1
l2: GROUP ATOMS=8                # Ligand atom 2
l3: GROUP ATOMS=7                # Ligand atom 3
l4: GROUP ATOMS=10               # Ligand atom 4

# System treatment
WHOLEMOLECULES ENTITY0=HOST
FIT_TO_TEMPLATE STRIDE=1 REFERENCE=host_template.pdb TYPE=OPTIMAL
lig: CENTER ATOMS=LIGC

# Virtual atoms defining funnel axis
v1: FIXEDATOM AT=0.000,0.000,0.000
v2: FIXEDATOM AT=0.000,0.000,0.250
v3: FIXEDATOM AT=0.000,0.000,0.500
v4: FIXEDATOM AT=0.000,0.000,0.750
v5: FIXEDATOM AT=0.000,0.000,1.000
v6: FIXEDATOM AT=0.000,0.000,-0.250
v7: FIXEDATOM AT=0.000,0.000,-0.500
v8: FIXEDATOM AT=0.000,0.000,-0.750
v9: FIXEDATOM AT=0.000,0.000,-1.000

# Basic position variables
cyl: DISTANCE ATOMS=v1,lig COMPONENTS
radius: MATHEVAL ARG=cyl.x,cyl.y FUNC=sqrt(x*x+y*y) PERIODIC=NO

# Ligand-water coordination
L1: COORDINATION GROUPA=l1 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
L2: COORDINATION GROUPA=l2 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
L3: COORDINATION GROUPA=l3 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
L4: COORDINATION GROUPA=l4 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20

# Virtual atom-water coordination
V1: COORDINATION GROUPA=v1 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V2: COORDINATION GROUPA=v2 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V3: COORDINATION GROUPA=v3 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V4: COORDINATION GROUPA=v4 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V5: COORDINATION GROUPA=v5 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V6: COORDINATION GROUPA=v6 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V7: COORDINATION GROUPA=v7 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V8: COORDINATION GROUPA=v8 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V9: COORDINATION GROUPA=v9 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
# Double funnel for symmetric hosts (e.g., CB8)
funnel: MATHEVAL ARG=radius,cyl.z VAR=r,z FUNC=((r+1.5*(-2.0+z))*step(-z+0.5)+(r-0.2)*step(z-0.5))*step(z-0.0)+((r-1.5*(2.0+z))*step(z+0.5)+(r-0.2)*step(-z-0.5))*step(-z+0.0) PERIODIC=NO

# Walls
UPPER_WALLS AT=0 ARG=funnel KAPPA=2000.0 LABEL=funnelwall
UPPER_WALLS AT=1.5 ARG=cyl.z KAPPA=2000.0 LABEL=upper_wall
LOWER_WALLS AT=-1.5 ARG=cyl.z KAPPA=2000.0 LABEL=lower_wall

# Angle and energy
ang: ANGLE ATOMS=v3,v5,1,10
cosang: MATHEVAL ARG=ang FUNC=cos(x) PERIODIC=NO
ene: ENERGY

OPES_METAD_EXPLORE ...
    LABEL=opes
    ARG=cyl.z,cosang
    SIGMA=0.1,0.2    # Use values from unbiased run
    FILE=Kernels.data
    STATE_RFILE=compressed.Kernels.data
    STATE_WFILE=compressed.Kernels.data
    PACE=10000
    BARRIER=100
... OPES_METAD_EXPLORE
OPES_METAD_EXPLORE ...
    LABEL=opese1
    ARG=L4
    SIGMA=0.3    # From unbiased run
    FILE=Kernels1.data
    STATE_RFILE=compressed_Kernels1.data
    STATE_WFILE=compressed.Kernels1.data
    PACE=20000
    BARRIER=3
... OPES_METAD_EXPLORE
OPES_METAD_EXPLORE ...
    LABEL=opese2
    ARG=V1
    SIGMA=0.3   # From unbiased run
    FILE=Kernels2.data
    STATE_RFILE=compressed_Kernels2.data
    STATE_WFILE=compressed.Kernels2.data
    PACE=20000
    BARRIER=3
... OPES_METAD_EXPLORE
#ENDHIDDEN
OPES_METAD_EXPLORE ...
    LABEL=opese3
    ARG=L1
    SIGMA=0.3   # From unbiased run
    FILE=Kernels3.data
    STATE_RFILE=compressed_Kernels3.data
    STATE_WFILE=compressed.Kernels3.data
    PACE=20000
    BARRIER=3
... OPES_METAD_EXPLORE

# Update print statement

PRINT ARG=opes.bias,cyl.z,radius,funnelwall.bias,upper_wall.bias,lower_wall.bias,ene,cosang,L1,L2,L3,L4,V1,V2,V3,V4,V5,V6,V7,V8,V9,opese1.bias,opese2.bias,opese3.bias STRIDE=50 FILE=COLVAR FMT=%8.4f
```
### Replica 4
Add ECV_MULTITHERMAL to previous replica:
```plumed
#HIDDEN
# Host, ligand and water groups
HOST: GROUP ATOMS=28-171         # Host atoms
LIGC: GROUP ATOMS=1-11           # Ligand heavy atoms
WO: GROUP ATOMS=172-4234:3       # Water oxygen atoms

# Selected ligand atoms for orientation and water coordination
l1: GROUP ATOMS=1                # Ligand atom 1
l2: GROUP ATOMS=8                # Ligand atom 2
l3: GROUP ATOMS=7                # Ligand atom 3
l4: GROUP ATOMS=10               # Ligand atom 4

# System treatment
WHOLEMOLECULES ENTITY0=HOST
FIT_TO_TEMPLATE STRIDE=1 REFERENCE=host_template.pdb TYPE=OPTIMAL
lig: CENTER ATOMS=LIGC

# Virtual atoms defining funnel axis
v1: FIXEDATOM AT=0.000,0.000,0.000
v2: FIXEDATOM AT=0.000,0.000,0.250
v3: FIXEDATOM AT=0.000,0.000,0.500
v4: FIXEDATOM AT=0.000,0.000,0.750
v5: FIXEDATOM AT=0.000,0.000,1.000
v6: FIXEDATOM AT=0.000,0.000,-0.250
v7: FIXEDATOM AT=0.000,0.000,-0.500
v8: FIXEDATOM AT=0.000,0.000,-0.750
v9: FIXEDATOM AT=0.000,0.000,-1.000

# Basic position variables
cyl: DISTANCE ATOMS=v1,lig COMPONENTS
radius: MATHEVAL ARG=cyl.x,cyl.y FUNC=sqrt(x*x+y*y) PERIODIC=NO

# Ligand-water coordination
L1: COORDINATION GROUPA=l1 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
L2: COORDINATION GROUPA=l2 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
L3: COORDINATION GROUPA=l3 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
L4: COORDINATION GROUPA=l4 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20

# Virtual atom-water coordination
V1: COORDINATION GROUPA=v1 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V2: COORDINATION GROUPA=v2 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V3: COORDINATION GROUPA=v3 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V4: COORDINATION GROUPA=v4 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V5: COORDINATION GROUPA=v5 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V6: COORDINATION GROUPA=v6 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V7: COORDINATION GROUPA=v7 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V8: COORDINATION GROUPA=v8 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V9: COORDINATION GROUPA=v9 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
# Double funnel for symmetric hosts (e.g., CB8)
funnel: MATHEVAL ARG=radius,cyl.z VAR=r,z FUNC=((r+1.5*(-2.0+z))*step(-z+0.5)+(r-0.2)*step(z-0.5))*step(z-0.0)+((r-1.5*(2.0+z))*step(z+0.5)+(r-0.2)*step(-z-0.5))*step(-z+0.0) PERIODIC=NO

# Walls
UPPER_WALLS AT=0 ARG=funnel KAPPA=2000.0 LABEL=funnelwall
UPPER_WALLS AT=1.5 ARG=cyl.z KAPPA=2000.0 LABEL=upper_wall
LOWER_WALLS AT=-1.5 ARG=cyl.z KAPPA=2000.0 LABEL=lower_wall

# Angle and energy
ang: ANGLE ATOMS=v3,v5,1,10
cosang: MATHEVAL ARG=ang FUNC=cos(x) PERIODIC=NO
ene: ENERGY
#ENDHIDDEN

# --- (4) OPES  ---
ecv: ECV_MULTITHERMAL ARG=ene TEMP_MAX=310
opesX: OPES_EXPANDED ARG=ecv.* FILE=DeltaFs.data PACE=100

#HIDDEN
OPES_METAD_EXPLORE ...
    LABEL=opes
    ARG=cyl.z,cosang
    SIGMA=0.1,0.2    # Use values from unbiased run
    FILE=Kernels.data
    STATE_RFILE=compressed.Kernels.data
    STATE_WFILE=compressed.Kernels.data
    PACE=10000
    BARRIER=100
... OPES_METAD_EXPLORE
OPES_METAD_EXPLORE ...
    LABEL=opese1
    ARG=L4
    SIGMA=0.3    # From unbiased run
    FILE=Kernels1.data
    STATE_RFILE=compressed_Kernels1.data
    STATE_WFILE=compressed.Kernels1.data
    PACE=20000
    BARRIER=3
... OPES_METAD_EXPLORE
OPES_METAD_EXPLORE ...
    LABEL=opese2
    ARG=V1
    SIGMA=0.3   # From unbiased run
    FILE=Kernels2.data
    STATE_RFILE=compressed_Kernels2.data
    STATE_WFILE=compressed.Kernels2.data
    PACE=20000
    BARRIER=3
... OPES_METAD_EXPLORE
OPES_METAD_EXPLORE ...
    LABEL=opese3
    ARG=L1
    SIGMA=0.3   # From unbiased run
    FILE=Kernels3.data
    STATE_RFILE=compressed_Kernels3.data
    STATE_WFILE=compressed.Kernels3.data
    PACE=20000
    BARRIER=3
... OPES_METAD_EXPLORE
#ENDHIDDEN

# Add CV to previous replica:
OPES_METAD_EXPLORE ...
   LABEL=opese4
   ARG=V3
   SIGMA=0.250   # From unbiased run
   FILE=Kernels4.data
   STATE_RFILE=compressed_Kernels4.data
   STATE_WFILE=compressed.Kernels4.data
   PACE=20000
   BARRIER=3
... OPES_METAD_EXPLORE

# Update print statement
PRINT ARG=opes.bias,cyl.z,radius,funnelwall.bias,upper_wall.bias,lower_wall.bias,ene,cosang,L1,L2,L3,L4,V1,V2,V3,V4,V5,V6,V7,V8,V9,opesX.bias,opese1.bias,opese2.bias,opese3.bias,opese4.bias STRIDE=50 FILE=COLVAR FMT=%8.4f
```
### Replica 5
Increase temperature of ECV_MULTITHERMAL:
```plumed
#HIDDEN
# Host, ligand and water groups
HOST: GROUP ATOMS=28-171         # Host atoms
LIGC: GROUP ATOMS=1-11           # Ligand heavy atoms
WO: GROUP ATOMS=172-4234:3       # Water oxygen atoms

# Selected ligand atoms for orientation and water coordination
l1: GROUP ATOMS=1                # Ligand atom 1
l2: GROUP ATOMS=8                # Ligand atom 2
l3: GROUP ATOMS=7                # Ligand atom 3
l4: GROUP ATOMS=10               # Ligand atom 4

# System treatment
WHOLEMOLECULES ENTITY0=HOST
FIT_TO_TEMPLATE STRIDE=1 REFERENCE=host_template.pdb TYPE=OPTIMAL
lig: CENTER ATOMS=LIGC

# Virtual atoms defining funnel axis
v1: FIXEDATOM AT=0.000,0.000,0.000
v2: FIXEDATOM AT=0.000,0.000,0.250
v3: FIXEDATOM AT=0.000,0.000,0.500
v4: FIXEDATOM AT=0.000,0.000,0.750
v5: FIXEDATOM AT=0.000,0.000,1.000
v6: FIXEDATOM AT=0.000,0.000,-0.250
v7: FIXEDATOM AT=0.000,0.000,-0.500
v8: FIXEDATOM AT=0.000,0.000,-0.750
v9: FIXEDATOM AT=0.000,0.000,-1.000

# Basic position variables
cyl: DISTANCE ATOMS=v1,lig COMPONENTS
radius: MATHEVAL ARG=cyl.x,cyl.y FUNC=sqrt(x*x+y*y) PERIODIC=NO

# Ligand-water coordination
L1: COORDINATION GROUPA=l1 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
L2: COORDINATION GROUPA=l2 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
L3: COORDINATION GROUPA=l3 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
L4: COORDINATION GROUPA=l4 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20

# Virtual atom-water coordination
V1: COORDINATION GROUPA=v1 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V2: COORDINATION GROUPA=v2 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V3: COORDINATION GROUPA=v3 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V4: COORDINATION GROUPA=v4 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V5: COORDINATION GROUPA=v5 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V6: COORDINATION GROUPA=v6 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V7: COORDINATION GROUPA=v7 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V8: COORDINATION GROUPA=v8 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V9: COORDINATION GROUPA=v9 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
# Double funnel for symmetric hosts (e.g., CB8)
funnel: MATHEVAL ARG=radius,cyl.z VAR=r,z FUNC=((r+1.5*(-2.0+z))*step(-z+0.5)+(r-0.2)*step(z-0.5))*step(z-0.0)+((r-1.5*(2.0+z))*step(z+0.5)+(r-0.2)*step(-z-0.5))*step(-z+0.0) PERIODIC=NO

# Walls
UPPER_WALLS AT=0 ARG=funnel KAPPA=2000.0 LABEL=funnelwall
UPPER_WALLS AT=1.5 ARG=cyl.z KAPPA=2000.0 LABEL=upper_wall
LOWER_WALLS AT=-1.5 ARG=cyl.z KAPPA=2000.0 LABEL=lower_wall

# Angle and energy
ang: ANGLE ATOMS=v3,v5,1,10
cosang: MATHEVAL ARG=ang FUNC=cos(x) PERIODIC=NO
ene: ENERGY
#ENDHIDDEN

# --- (4) OPES  ---
ecv: ECV_MULTITHERMAL ARG=ene TEMP_MAX=330
opesX: OPES_EXPANDED ARG=ecv.* FILE=DeltaFs.data PACE=100

#HIDDEN
OPES_METAD_EXPLORE ...
    LABEL=opes
    ARG=cyl.z,cosang
    SIGMA=0.1,0.2    # Use values from unbiased run
    FILE=Kernels.data
    STATE_RFILE=compressed.Kernels.data
    STATE_WFILE=compressed.Kernels.data
    PACE=10000
    BARRIER=100
... OPES_METAD_EXPLORE
OPES_METAD_EXPLORE ...
    LABEL=opese1
    ARG=L4
    SIGMA=0.3    # From unbiased run
    FILE=Kernels1.data
    STATE_RFILE=compressed_Kernels1.data
    STATE_WFILE=compressed.Kernels1.data
    PACE=20000
    BARRIER=3
... OPES_METAD_EXPLORE
OPES_METAD_EXPLORE ...
    LABEL=opese2
    ARG=V1
    SIGMA=0.3   # From unbiased run
    FILE=Kernels2.data
    STATE_RFILE=compressed_Kernels2.data
    STATE_WFILE=compressed.Kernels2.data
    PACE=20000
    BARRIER=3
... OPES_METAD_EXPLORE
OPES_METAD_EXPLORE ...
    LABEL=opese3
    ARG=L1
    SIGMA=0.3   # From unbiased run
    FILE=Kernels3.data
    STATE_RFILE=compressed_Kernels3.data
    STATE_WFILE=compressed.Kernels3.data
    PACE=20000
    BARRIER=3
... OPES_METAD_EXPLORE
# Add CV to previous replica:
OPES_METAD_EXPLORE ...
   LABEL=opese4
   ARG=V3
   SIGMA=0.250   # From unbiased run
   FILE=Kernels4.data
   STATE_RFILE=compressed_Kernels4.data
   STATE_WFILE=compressed.Kernels4.data
   PACE=20000
   BARRIER=3
... OPES_METAD_EXPLORE
#ENDHIDDEN

# Add CV to previous replica:
OPES_METAD_EXPLORE ...
   LABEL=opese5
   ARG=V7
   SIGMA=0.440   # From unbiased run
   FILE=Kernels5.data
   STATE_RFILE=compressed_Kernels5.data
   STATE_WFILE=compressed.Kernels5.data
   PACE=20000
   BARRIER=3
... OPES_METAD_EXPLORE

# Update print statement
PRINT ARG=opes.bias,cyl.z,radius,funnelwall.bias,upper_wall.bias,lower_wall.bias,ene,cosang,L1,L2,L3,L4,V1,V2,V3,V4,V5,V6,V7,V8,V9,opesX.bias,opese1.bias,opese2.bias,opese3.bias,opese4.bias,opese5.bias STRIDE=50 FILE=COLVAR FMT=%8.4f
```
### Replica 6
Increase temperature:
```plumed
#HIDDEN
# Host, ligand and water groups
HOST: GROUP ATOMS=28-171         # Host atoms
LIGC: GROUP ATOMS=1-11           # Ligand heavy atoms
WO: GROUP ATOMS=172-4234:3       # Water oxygen atoms

# Selected ligand atoms for orientation and water coordination
l1: GROUP ATOMS=1                # Ligand atom 1
l2: GROUP ATOMS=8                # Ligand atom 2
l3: GROUP ATOMS=7                # Ligand atom 3
l4: GROUP ATOMS=10               # Ligand atom 4

# System treatment
WHOLEMOLECULES ENTITY0=HOST
FIT_TO_TEMPLATE STRIDE=1 REFERENCE=host_template.pdb TYPE=OPTIMAL
lig: CENTER ATOMS=LIGC

# Virtual atoms defining funnel axis
v1: FIXEDATOM AT=0.000,0.000,0.000
v2: FIXEDATOM AT=0.000,0.000,0.250
v3: FIXEDATOM AT=0.000,0.000,0.500
v4: FIXEDATOM AT=0.000,0.000,0.750
v5: FIXEDATOM AT=0.000,0.000,1.000
v6: FIXEDATOM AT=0.000,0.000,-0.250
v7: FIXEDATOM AT=0.000,0.000,-0.500
v8: FIXEDATOM AT=0.000,0.000,-0.750
v9: FIXEDATOM AT=0.000,0.000,-1.000

# Basic position variables
cyl: DISTANCE ATOMS=v1,lig COMPONENTS
radius: MATHEVAL ARG=cyl.x,cyl.y FUNC=sqrt(x*x+y*y) PERIODIC=NO

# Ligand-water coordination
L1: COORDINATION GROUPA=l1 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
L2: COORDINATION GROUPA=l2 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
L3: COORDINATION GROUPA=l3 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
L4: COORDINATION GROUPA=l4 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20

# Virtual atom-water coordination
V1: COORDINATION GROUPA=v1 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V2: COORDINATION GROUPA=v2 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V3: COORDINATION GROUPA=v3 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V4: COORDINATION GROUPA=v4 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V5: COORDINATION GROUPA=v5 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V6: COORDINATION GROUPA=v6 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V7: COORDINATION GROUPA=v7 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V8: COORDINATION GROUPA=v8 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V9: COORDINATION GROUPA=v9 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
# Double funnel for symmetric hosts (e.g., CB8)
funnel: MATHEVAL ARG=radius,cyl.z VAR=r,z FUNC=((r+1.5*(-2.0+z))*step(-z+0.5)+(r-0.2)*step(z-0.5))*step(z-0.0)+((r-1.5*(2.0+z))*step(z+0.5)+(r-0.2)*step(-z-0.5))*step(-z+0.0) PERIODIC=NO

# Walls
UPPER_WALLS AT=0 ARG=funnel KAPPA=2000.0 LABEL=funnelwall
UPPER_WALLS AT=1.5 ARG=cyl.z KAPPA=2000.0 LABEL=upper_wall
LOWER_WALLS AT=-1.5 ARG=cyl.z KAPPA=2000.0 LABEL=lower_wall

# Angle and energy
ang: ANGLE ATOMS=v3,v5,1,10
cosang: MATHEVAL ARG=ang FUNC=cos(x) PERIODIC=NO
ene: ENERGY
#ENDHIDDEN

# --- (4) OPES  ---
ecv: ECV_MULTITHERMAL ARG=ene TEMP_MAX=350
opesX: OPES_EXPANDED ARG=ecv.* FILE=DeltaFs.data PACE=100

#HIDDEN
OPES_METAD_EXPLORE ...
    LABEL=opes
    ARG=cyl.z,cosang
    SIGMA=0.1,0.2    # Use values from unbiased run
    FILE=Kernels.data
    STATE_RFILE=compressed.Kernels.data
    STATE_WFILE=compressed.Kernels.data
    PACE=10000
    BARRIER=100
... OPES_METAD_EXPLORE
OPES_METAD_EXPLORE ...
    LABEL=opese1
    ARG=L4
    SIGMA=0.3    # From unbiased run
    FILE=Kernels1.data
    STATE_RFILE=compressed_Kernels1.data
    STATE_WFILE=compressed.Kernels1.data
    PACE=20000
    BARRIER=3
... OPES_METAD_EXPLORE
OPES_METAD_EXPLORE ...
    LABEL=opese2
    ARG=V1
    SIGMA=0.3   # From unbiased run
    FILE=Kernels2.data
    STATE_RFILE=compressed_Kernels2.data
    STATE_WFILE=compressed.Kernels2.data
    PACE=20000
    BARRIER=3
... OPES_METAD_EXPLORE
OPES_METAD_EXPLORE ...
    LABEL=opese3
    ARG=L1
    SIGMA=0.3   # From unbiased run
    FILE=Kernels3.data
    STATE_RFILE=compressed_Kernels3.data
    STATE_WFILE=compressed.Kernels3.data
    PACE=20000
    BARRIER=3
... OPES_METAD_EXPLORE
# Add CV to previous replica:
OPES_METAD_EXPLORE ...
   LABEL=opese4
   ARG=V3
   SIGMA=0.250   # From unbiased run
   FILE=Kernels4.data
   STATE_RFILE=compressed_Kernels4.data
   STATE_WFILE=compressed.Kernels4.data
   PACE=20000
   BARRIER=3
... OPES_METAD_EXPLORE
# Add CV to previous replica:
OPES_METAD_EXPLORE ...
   LABEL=opese5
   ARG=V7
   SIGMA=0.440   # From unbiased run
   FILE=Kernels5.data
   STATE_RFILE=compressed_Kernels5.data
   STATE_WFILE=compressed.Kernels5.data
   PACE=20000
   BARRIER=3
... OPES_METAD_EXPLORE
#ENDHIDDEN

# Add CV to previous replica:
OPES_METAD_EXPLORE ...
   LABEL=opese6
   ARG=V5
   SIGMA=0.280   # From unbiased run
   FILE=Kernels6.data
   STATE_RFILE=compressed_Kernels6.data
   STATE_WFILE=compressed.Kernels6.data
   PACE=20000
   BARRIER=3
... OPES_METAD_EXPLORE

# Update print statement
PRINT ARG=opes.bias,cyl.z,radius,funnelwall.bias,upper_wall.bias,lower_wall.bias,ene,cosang,L1,L2,L3,L4,V1,V2,V3,V4,V5,V6,V7,V8,V9,opesX.bias,opese1.bias,opese2.bias,opese3.bias,opese4.bias,opese5.bias,opese6.bias STRIDE=50 FILE=COLVAR FMT=%8.4f
```
### Replica 7
Increase temperature:
```plumed
#HIDDEN
# Host, ligand and water groups
HOST: GROUP ATOMS=28-171         # Host atoms
LIGC: GROUP ATOMS=1-11           # Ligand heavy atoms
WO: GROUP ATOMS=172-4234:3       # Water oxygen atoms

# Selected ligand atoms for orientation and water coordination
l1: GROUP ATOMS=1                # Ligand atom 1
l2: GROUP ATOMS=8                # Ligand atom 2
l3: GROUP ATOMS=7                # Ligand atom 3
l4: GROUP ATOMS=10               # Ligand atom 4

# System treatment
WHOLEMOLECULES ENTITY0=HOST
FIT_TO_TEMPLATE STRIDE=1 REFERENCE=host_template.pdb TYPE=OPTIMAL
lig: CENTER ATOMS=LIGC

# Virtual atoms defining funnel axis
v1: FIXEDATOM AT=0.000,0.000,0.000
v2: FIXEDATOM AT=0.000,0.000,0.250
v3: FIXEDATOM AT=0.000,0.000,0.500
v4: FIXEDATOM AT=0.000,0.000,0.750
v5: FIXEDATOM AT=0.000,0.000,1.000
v6: FIXEDATOM AT=0.000,0.000,-0.250
v7: FIXEDATOM AT=0.000,0.000,-0.500
v8: FIXEDATOM AT=0.000,0.000,-0.750
v9: FIXEDATOM AT=0.000,0.000,-1.000

# Basic position variables
cyl: DISTANCE ATOMS=v1,lig COMPONENTS
radius: MATHEVAL ARG=cyl.x,cyl.y FUNC=sqrt(x*x+y*y) PERIODIC=NO

# Ligand-water coordination
L1: COORDINATION GROUPA=l1 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
L2: COORDINATION GROUPA=l2 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
L3: COORDINATION GROUPA=l3 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
L4: COORDINATION GROUPA=l4 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=6 MM=10 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20

# Virtual atom-water coordination
V1: COORDINATION GROUPA=v1 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V2: COORDINATION GROUPA=v2 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V3: COORDINATION GROUPA=v3 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V4: COORDINATION GROUPA=v4 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V5: COORDINATION GROUPA=v5 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V6: COORDINATION GROUPA=v6 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V7: COORDINATION GROUPA=v7 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V8: COORDINATION GROUPA=v8 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
V9: COORDINATION GROUPA=v9 GROUPB=WO SWITCH={RATIONAL D_0=0.0 R_0=0.25 NN=2 MM=6 D_MAX=0.8} NLIST NL_CUTOFF=1.5 NL_STRIDE=20
# Double funnel for symmetric hosts (e.g., CB8)
funnel: MATHEVAL ARG=radius,cyl.z VAR=r,z FUNC=((r+1.5*(-2.0+z))*step(-z+0.5)+(r-0.2)*step(z-0.5))*step(z-0.0)+((r-1.5*(2.0+z))*step(z+0.5)+(r-0.2)*step(-z-0.5))*step(-z+0.0) PERIODIC=NO

# Walls
UPPER_WALLS AT=0 ARG=funnel KAPPA=2000.0 LABEL=funnelwall
UPPER_WALLS AT=1.5 ARG=cyl.z KAPPA=2000.0 LABEL=upper_wall
LOWER_WALLS AT=-1.5 ARG=cyl.z KAPPA=2000.0 LABEL=lower_wall

# Angle and energy
ang: ANGLE ATOMS=v3,v5,1,10
cosang: MATHEVAL ARG=ang FUNC=cos(x) PERIODIC=NO
ene: ENERGY
#ENDHIDDEN

# --- (4) OPES  ---
ecv: ECV_MULTITHERMAL ARG=ene TEMP_MAX=370
opesX: OPES_EXPANDED ARG=ecv.* FILE=DeltaFs.data PACE=100

#HIDDEN
OPES_METAD_EXPLORE ...
    LABEL=opes
    ARG=cyl.z,cosang
    SIGMA=0.1,0.2    # Use values from unbiased run
    FILE=Kernels.data
    STATE_RFILE=compressed.Kernels.data
    STATE_WFILE=compressed.Kernels.data
    PACE=10000
    BARRIER=100
... OPES_METAD_EXPLORE
OPES_METAD_EXPLORE ...
    LABEL=opese1
    ARG=L4
    SIGMA=0.3    # From unbiased run
    FILE=Kernels1.data
    STATE_RFILE=compressed_Kernels1.data
    STATE_WFILE=compressed.Kernels1.data
    PACE=20000
    BARRIER=3
... OPES_METAD_EXPLORE
OPES_METAD_EXPLORE ...
    LABEL=opese2
    ARG=V1
    SIGMA=0.3   # From unbiased run
    FILE=Kernels2.data
    STATE_RFILE=compressed_Kernels2.data
    STATE_WFILE=compressed.Kernels2.data
    PACE=20000
    BARRIER=3
... OPES_METAD_EXPLORE
OPES_METAD_EXPLORE ...
    LABEL=opese3
    ARG=L1
    SIGMA=0.3   # From unbiased run
    FILE=Kernels3.data
    STATE_RFILE=compressed_Kernels3.data
    STATE_WFILE=compressed.Kernels3.data
    PACE=20000
    BARRIER=3
... OPES_METAD_EXPLORE
# Add CV to previous replica:
OPES_METAD_EXPLORE ...
   LABEL=opese4
   ARG=V3
   SIGMA=0.250   # From unbiased run
   FILE=Kernels4.data
   STATE_RFILE=compressed_Kernels4.data
   STATE_WFILE=compressed.Kernels4.data
   PACE=20000
   BARRIER=3
... OPES_METAD_EXPLORE
# Add CV to previous replica:
OPES_METAD_EXPLORE ...
   LABEL=opese5
   ARG=V7
   SIGMA=0.440   # From unbiased run
   FILE=Kernels5.data
   STATE_RFILE=compressed_Kernels5.data
   STATE_WFILE=compressed.Kernels5.data
   PACE=20000
   BARRIER=3
... OPES_METAD_EXPLORE
# Add CV to previous replica:
OPES_METAD_EXPLORE ...
   LABEL=opese6
   ARG=V5
   SIGMA=0.280   # From unbiased run
   FILE=Kernels6.data
   STATE_RFILE=compressed_Kernels6.data
   STATE_WFILE=compressed.Kernels6.data
   PACE=20000
   BARRIER=3
... OPES_METAD_EXPLORE
#ENDHIDDEN

# Add CV to previous replica:
OPES_METAD_EXPLORE ...
   LABEL=opese7
   ARG=V9
   SIGMA=0.630   # From unbiased run
   FILE=Kernels7.data
   STATE_RFILE=compressed_Kernels7.data
   STATE_WFILE=compressed.Kernels7.data
   PACE=20000
   BARRIER=3
... OPES_METAD_EXPLORE

# Update print statement
PRINT ARG=opes.bias,cyl.z,radius,funnelwall.bias,upper_wall.bias,lower_wall.bias,ene,cosang,L1,L2,L3,L4,V1,V2,V3,V4,V5,V6,V7,V8,V9,opesX.bias,opese1.bias,opese2.bias,opese3.bias,opese4.bias,opese5.bias,opese6.bias,opese7.bias STRIDE=50 FILE=COLVAR FMT=%8.4f
```
## 5. Running OneOPES Simulation
Create 8 directories, named from 0 to 7. Copy the tpr, topology and gro files, the host_template.pdb and the corresponding plumed.dat file to each folder. Then run OneOpes multireplica with the displayed mpirun command. 
```bash
# Create directories
mkdir {0..7}
cp topol.top {0..7}/
cp npt.gro {0..7}/
cp topol.tpr {0..7}/
cp host_template.pdb {0..7}/
cp plumed_X.dat X/plumed.dat  # For each replica X

# Run simulation
mpirun -n 8 gmx_mpi mdrun -deffnm topol -multidir 0 1 2 3 4 5 6 7 -replex 1000 -hrex -plumed plumed.dat -nsteps 100000000 -notunepme
```

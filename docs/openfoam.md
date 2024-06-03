# OpenFOAM

[OpenFOAM](http://www.openfoam.org) (Open Field Operation and Manipulation) is an open-source computational fluid dynamics (CFD) software used for simulating fluid flow, turbulence, heat transfer, and reactions in gases, liquids, and solids. It is designed for flexibility, extensibility, and high performance, offering a wide range of features for modeling and solving complex fluid dynamics problems.

## Key Components of OpenFOAM

1. **Mesh Generation:**
   - **BlockMesh:** A basic mesh generator using simple geometrical shapes.
   - **SnappyHexMesh:** A more advanced tool for generating hexahedral-dominant meshes that conform to complex geometries.
   - **TopoSet:** A utility for creating and manipulating sets of mesh entities such as cells, faces, and points.

2. **Solvers:**
   - OpenFOAM includes a variety of solvers for different types of simulations, such as incompressible and compressible flow, multiphase flow, heat transfer, and chemical reactions.
   - Users can also develop custom solvers by modifying the existing ones or creating new ones from scratch using OpenFOAM's extensive C++ library.

3. **Pre-Processing:**
   - Preparation of input files, such as setting boundary conditions, initial conditions, and physical properties.
   - Configuration of dictionaries for different simulation aspects, such as mesh generation (`blockMeshDict`, `snappyHexMeshDict`) and solver settings (`controlDict`, `fvSchemes`, `fvSolution`).

4. **Simulation Execution:**
   - Running simulations by executing the appropriate solver.
   - Parallel computation capabilities allow for efficient simulations on multi-core processors and high-performance computing clusters.

5. **Post-Processing:**
   - Tools like ParaView for visualizing simulation results.
   - Utilities for data manipulation, such as sampling and processing fields, generating reports, and plotting results.

## How OpenFOAM Works

1. **Problem Setup:**
   - Define the problem domain, including geometry and mesh generation.
   - Set up initial and boundary conditions.
   - Configure the physical properties of the fluid or solid.

2. **Mesh Generation:**
   - Use BlockMesh or SnappyHexMesh to create a computational mesh that discretizes the problem domain.
   - Refine and adapt the mesh to accurately capture the geometry and important flow features.

3. **Define Simulation Parameters:**
   - Configure solver settings and numerical schemes.
   - Specify control parameters such as time stepping, convergence criteria, and output intervals.

4. **Run the Simulation:**
   - Execute the solver to simulate the physical phenomena over the defined time period or until convergence.
   - Monitor the simulation progress and adjust parameters if necessary to ensure accuracy and stability.

5. **Post-Process Results:**
   - Visualize the flow fields, temperature distribution, or other relevant quantities using visualization tools.
   - Analyze the results to extract meaningful insights, such as pressure drops, velocity profiles, or heat fluxes.

## Advantages of OpenFOAM

- **Open Source:** Free to use and modify, with a large user community contributing to its development and support.
- **Extensibility:** Highly customizable, allowing users to create new solvers, utilities, and libraries.
- **Flexibility:** Supports a wide range of physical models and can be applied to various engineering problems.
- **Parallel Processing:** Efficient parallel computation capabilities for large-scale simulations.

## Challenges

- **Complexity:** Steep learning curve due to the extensive use of C++ and the need for familiarity with CFD concepts.
- **Documentation:** While comprehensive, the documentation can sometimes be fragmented or outdated.


# Installing OpenFOAM

To install check (https://openfoam.org/download/) and follow instructions.


## Prerequisites
- **Operating System:** Linux (e.g., Ubuntu), macOS, or Windows (via Windows Subsystem for Linux).
- **Dependencies:** GCC, CMake, Python, and additional libraries such as Scotch and ParaView for visualization.

## Installation Steps

### On Ubuntu
1. **Add the OpenFOAM repository:**
```bash
   sudo sh -c "wget -O - https://dl.openfoam.org/gpg.key | apt-key add -"
   sudo add-apt-repository http://dl.openfoam.org/ubuntu
```

2. **Update the package list:**
```bash
    sudo apt-get update
```
3. **Install OpenFOAM:**
```bash
    sudo apt-get install openfoam10
```
4. **Set up the environment:**
Add the following lines to your `~/.bashrc` file:
```bash
    source /opt/openfoam10/etc/bashrc
```
Then run
```bash
    source ~/.bashrc
```

### On MacOS

1. **Installing Docker for Mac**

    Follow Instructions in https://docs.docker.com/get-docker/

2. **Download openfoam10-macos**

```bash
    sudo curl --create-dirs -o /usr/local/bin/openfoam10-macos http://dl.openfoam.org/docker/openfoam10-macos
    sudo chmod 755 /usr/local/bin/openfoam10-macos
```

3. **Launching openfoam10-macos**

    The Docker container mounts the user’s file system so that case files are stored permanently. The container mounts the directory from where openfoam10-macos is launched by default

```bash
    cd $HOME/openfoam
    openfoam10-macos
```
### On Windows

The packaged distributions of OpenFOAM for Ubuntu can now be installed directly on Microsoft Windows 10 using Windows Subsystem for Linux (WSL).


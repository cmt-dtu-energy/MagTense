Demag field - FMM
=================

Calculating the demagnetization (stray) field is typically the most computationally demanding part of a micromagnetic simulation. While MagTense’s default analytical approach is exact, its $O(N^2)$ scaling can be prohibitive for very large systems. 

To address this, MagTense includes an implementation of the **Fast Multipole Method (FMM)**, which reduces the computational complexity to $O(N)$.

Implementation & Attribution
---------------------------
The FMM acceleration in MagTense is built upon the **FMM3D** library developed by the **FlatIron Institute**. 

* **Core Library:** `FlatIron Institute FMM3D <https://github.com/flatironinstitute/fmm3d>`_
* **Linking:** MagTense is linked with a specialized fork, `Ximtecs/FMM3D <https://github.com/Ximtecs/FMM3D>`_, which contains updated build configurations and makefiles to allow seamless integration with the MagTense Fortran core.

Compilation
-----------
To enable FMM support during the build process, you must explicitly include the FMM flag in your make command. This ensures the compiler links the necessary FMM3D libraries and enables the specialized Fortran modules.

.. code-block:: bash

    make USE_FMM3D=1

Optimization: Persistent Tree Structure
--------------------------------------
In micromagnetic simulations, the spatial distribution of cells (the mesh) is typically static throughout the simulation. To maximize efficiency, MagTense utilizes a **"magtense-local"** tree structure. 

This structure caches the octree setup between consecutive calls to the solver. By reusing the tree, the overhead of re-partitioning space is eliminated for every time step, significantly accelerating the total simulation time. 

.. note::
   The current implementation requires a fully grown tree. Therefore, **ifunif** must always be set to **1**.

Near-Field Evaluation & Neighbor Tensors
----------------------------------------
In standard FMM implementations, near-field interactions are usually handled by a direct "Point-to-Point" (P2P) evaluation. In MagTense, this standard P2P evaluation is **disabled**. 

Instead, MagTense utilizes its high-precision analytical demagnetization tensors for all near-field interactions:

1. **Neighbor Identification:** Based on **List 1** from the FMM tree creation, MagTense identifies "neighbor" pairs (cells within the same or adjacent leaf nodes).
2. **Sparse Neighbor Tensor:** A sparse tensor structure is created to map these neighbor interactions.
3. **Analytical Calculation:** The exact analytical demagnetization tensor is calculated for every neighbor pair.
4. **Sparse Matrix Storage:** These values are converted into a sparse matrix. 
   - On **CPU**, this is stored as an **Intel MKL sparse matrix**.
   - On **GPU**, it is stored persistently in global memory for **CUDA** evaluation in each timestep.

Input Variables
---------------

.. list-table::
   :widths: 25 10 65
   :header-rows: 1

   * - Variable
     - Type
     - Description
   * - **fmm_cells_per_node**
     - int
     - Maximum cells in a leaf node before splitting.
   * - **eps_fmm**
     - float
     - Controls the **exponential order** (plane-wave expansion).
   * - **fmm_nterms**
     - int
     - Sets the **multipole expansion order**. (0 = auto).
   * - **ifunif**
     - int
     - **Tree Type:** Must be **1** (Uniform tree).
   * - **nlmin**
     - int
     - Minimum level of the octree hierarchy.
   * - **nlmax**
     - int
     - Maximum level (dictates the depth of the uniform tree).
   * - **allow_fmm_short_circuit**
     - int
     - **1** to allow direct calculation for small problems; **0** to force FMM.
   * - **fmm_min_n**
     - int
     - Cell count threshold required to activate FMM.

Python Usage Example
--------------------

.. code-block:: python

    problem.fmm_cells_per_node = 10
    problem.eps_fmm = 1e-4          
    problem.fmm_nterms = 12         
    problem.ifunif = 1              
    problem.nlmin = 1
    problem.nlmax = 2
    problem.allow_fmm_short_circuit = 1
    problem.fmm_min_n = 20000        

    RunMicroMagSimulation(problem)
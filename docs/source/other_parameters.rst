Other parameters
========================================
::
  
    real,dimension(3) :: symmetryOps
    integer :: exploitSymmetry
    integer :: includeInIteration
    logical :: excludeFromSummation
    real,dimension(3) :: color

========================================
Symmetry Options
========================================
* 1 = symmetry
* -1 = for anti-symmetry 

| **symmetryOps**:math:`_1` for symmetry about xy plane.
| **symmetryOps**:math:`_2` for symmetry about xz plane.
| **symmetryOps**:math:`_3` for symmetry about yz plane.

For applying the given symmetry in the simulation, the parameter
**exploitSymmetry** has to be set to 1. Its default value is 0.

========================================
Include in Iteration
========================================
The parameter defines if the given tile is considered during the
iteration for finding the magnetization of each tile and
their effect among each other.
It has to be set to 0 for not including the current tile.
Its default value is set 1.

========================================
Exclude from summation
========================================
The parameter defines if the given tile should be included or 
not in the summation over all tiles for getting the applied H-field.   


========================================
Color
========================================
The color of the tile as rgb triplet.
It is used for displaying the model.

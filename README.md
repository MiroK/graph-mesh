# Meshing of graphs for FEniCS

## TODOs
-[ ] orientation of terminal nodes for in/outflow as determined by `TangentCurve`
-[ ] seperate graph generators (Alex and something more branching and loopy)
-[ ] what are `FEniCS_ii` extensions needed for [mixed Darcy](https://mox.polimi.it/reports-and-theses/publication-results/?id=632) or [Stokes](https://arxiv.org/abs/2111.12451)
-[ ] for presentations reduce number of colors used in `color_branches`

## Installation
For now run in bash `source setup.rc`

## Dependencies
- FEniCS stack
- graph computations require `networkx`
- some features depend `FEniCS_ii`
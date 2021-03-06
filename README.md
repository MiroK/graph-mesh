# Meshing of graphs for FEniCS

For usage see the [demo folder](https://github.com/MiroK/graph-mesh/blob/master/demo). For example
using the [BRAVA](http://cng.gmu.edu/brava/all_subjects.php?clear=1) database we can do something like this

  <p align="center">
    <img src="https://github.com/MiroK/graph-mesh/blob/master/doc/coloring.png">
  </p>

  <p align="center">
    <img src="https://github.com/MiroK/graph-mesh/blob/master/doc/tangent.png">
  </p>

## Installation
For now run in bash `source setup.rc`

## Dependencies
- FEniCS stack
- graph computations require `networkx`
- some features depend `FEniCS_ii`
- `embedding` uses [gmshnics](https://github.com/MiroK/gmshnics)

## TODOs
- [X] orientation of terminal nodes for in/outflow as determined by `TangentCurve`
- [X] seperate graph generators (Alex and something more branching and loopy)
- [ ] what are `FEniCS_ii` extensions needed for [mixed Darcy](https://mox.polimi.it/reports-and-theses/publication-results/?id=632) or [Stokes](https://arxiv.org/abs/2111.12451)
- [x] for presentations reduce number of colors used in `color_branches`
- [ ] align the bounding box in embedding
- [ ] allow for bounding volume as STL
- [ ] with boxes we connect the 1d to boundary (keep track of where)
- [ ] pass out radius and tagging info
# A Simple Event Tree for Studying Wire Cell Signals

## Instructions

### Build the Application
    cd scripts
    root -b -q -l loadClasses.C

- The simulated event is saved as a ROOT TTree object `CellTree`. The TTree is created with a LarSoft module [CellTree](CellTree).

- Example of how to read the event can be found at [scripts/example.C](scripts/example.C).

- Data files for test can be found at [data/filelist.txt](data/filelist.txt).

- Description of the data file version can be fount at [data/versions.md](data/versions.md).

- Geometry file for conversion can be fount at [geometry/ChannelWireGeometry.txt](geometry/ChannelWireGeometry.txt).

- The Geometry file is dumped from LarSoft using the `CTree::saveChannelWireGeometry()` method in [CTree/CTree_module.cc](CTree/CTree_module.cc)

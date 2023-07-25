# Changelog & Release Notes

## Upgrading

To upgrade to the latest version of `cfdtools` use `pip`:

```bash
pip install cfdtools --upgrade
```

You can determine your currently installed version using this command:

```bash
pip show cfdtools
```

## Versions

### [0.5.0](https://pypi.org/project/cfdtools/) (2023-xx-xx)

#### new

- Legendre polynomial extrapolation (development)
- DataSet and DataSetList classes
- read list of vtk files to DataSetList
- dump DataSetList as hdf5 file
- new `vtkbrief` command line tool

#### fix

- `ic3brief` no longer needs `--fmt IC3` for any file extension

### [0.4.2](https://pypi.org/project/cfdtools/) (2023-04-05)

#### new

- cgns import (hdf5 direct reader) (validated from icem and autogrid)
- extrusion of 2D mesh
- `cfdwrite_vtk` command line
- new option `--extrude n` to command lines
- new option `--scale x y z` to command lines
  
### [0.3.3](https://pypi.org/project/cfdtools/) (2023-03-13)

#### new

- internal structured 3D mesh `cfdtools.meshbase.simple.Cube`
- `cfdwritecube` command line
- support GMSH 2.x and 4.x as a reader
- `morph`function for meshbase object
- core: computation of face/node and face/cell connectivities from cell/node
- command `ic3probe_plotline`
- new common option `--outpath`
- export to vtk (not yet available as a command line)

#### changed

- new connectivities, boundary conditions marks

### [0.2.0](https://pypi.org/project/cfdtools/) (2022-03-22)

#### new

- command line `ic3brief` only section headers
- command line writers can remove date with `--remove-[cell/node]-data` options

#### changed

- ensure safe new names for writers to avoid deletion of read files

### [0.1.0](https://pypi.org/project/cfdtools/) (2022-03-11)

#### new

- command line `cfdwrite_ic3v2`, `cfdwrite_ic3v3`
- new `v3 IC3` writer

#### changed

- improve IC3 v2 and v3 with parameters and actual partition

### [0.0.2](https://pypi.org/project/cfdtools/) (2022-03-04)

#### new

- command line `cfdinfo`
- new `v2 IC3` writer

#### fixed

- v2 IC3 reader 

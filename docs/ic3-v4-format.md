
# IC3 V4 Format

Only for finite-volume computation.

The IC3 restart files use the HDF5 file format to manage, process, and store the restart data.
The data format is proprietary.

The mesh file is separated from the solution files when the mesh is not deformable.

HDF5 can store two primary types of objects: datasets and groups. A dataset is essentially a multidimensional array of data elements, and a group is a structure for organizing objects in an HDF5 file. The format supports user-defined attributes that can be used to add descriptive metadata to the file as a whole or any component data object.

IC3 group names start with a capital letter and have a capital letter for each new word, with no underscores except for groups related to specific labels.

IC3 dataset names are all lowercase.

## Mesh File Organization

The mesh filename should be `<mesh>.h5`.

### File Attributes

|   |   |   |
|---|---|---|
| ic3_format_version | 64-bit integer | optional  |
| hdf5_version | string | optional |
| cfdtools_version | string | optional |
| meshtype | string | optional |
| cv_count | 64-bit integer | required |
| fa_count | 64-bit integer | required |
| no_count | 64-bit integer | required |

### File Objects

|   |   |   |
|---|---|---|
| coordinates | HDF5 Dataset | required |
| Connectivities | HDF5 Group | required |
| Boundaries | HDF5 Group | required |
| NodeData | HDF5 Group | optional |

### Dataset `coordinates`

Point coordinates of the mesh.

Data type: 64-bit floating-point

Attributes:
|   |   |   |
|---|---|---|
| system | string | optional in ["cartesian", "cylindrical"] |

### Group `Connectivities`

|   |   |   |
|---|---|---|
| Face | HDF5 Group | required |
| Cell | HDF5 Group | optional |

The Face group stores the datasets for face-based connectivity.

|   |   |   |   |
|---|---|---|---|
| noofa_v | 32-bit integer | required | NGON_n (UGP_IO_NOOFA_I_AND_V, conn->noofa_v) |
| noofa_i | 32-bit integer | required | ElementStartOffset(NGON_n)  (kind of UGP_IO_NOOFA_I_AND_V, kind of conn->noofa_i) |
| cvofa | 32-bit integer | required | ParentElements (NGON_n) (UGP_IO_CVOFA, conn->cvofa) |

The Cell group stores the datasets for element-based connectivity.

|   |   |   |   |
|---|---|---|---|
| hexa8 | 64-bit integer | optional | (element-nodes connectivity) |

### Group `Boundaries`

|   |   |   |
|---|---|---|
| kinds | HDF5 Dataset, 32-bit integer | required |
| labels | HDF5 Dataset, string | required |
| offsets | HDF5 Dataset, 32-bit integer | required |

The kinds dataset stores the IC3 type of boundaries.

The offsets dataset stores the offset to reach the first face of the boundary in the global list of faces.

The labels dataset stores the name of the boundaries.

For each boundary name in the `labels` array, a HDF5 group with the exact same name may be present in the Boundaries group. If so, it contains a dataset to give the periodic transformation.

|   |   |   |
|---|---|---|
| < boundary name > | HDF5 Dataset, 64-bit floating-point | required |

### Group `NodeData`

|   |   |   |
|---|---|---|
| node_gb_index | HDF5 Dataset, 64-bit integer | optional |

## Solution File Organization

The solution filename should be `<restart>.<iteration>.h5`.

A solution file stores data that corresponds to a single physical time. The name should contain the iteration number of the computation corresponding to this time.

### File Attributes

|   |   |   |
|---|---|---|
| step | 64-bit integer | required |
| time | 64-bit floating-point | required |

The step attribute is the number of iterations already performed during previous runs.

The time attribute is the physical time of the stored data.

### File Objects

|   |   |   |
|---|---|---|
| rho | HDF5 Dataset | required |
| rhoe | HDF5 Dataset | required |
| rhou | HDF5 Dataset | required |
| VolumeStats | HDF5 Group | required |

### Datasets `rho`, `rhoe`, `rhou`

Conservative variables located at the control volume centers.

Data type: 64-bit floating-point

### Group `VolumeStats`

Thos group stores the time-averaged variables.

|   |   |   |
|---|---|---|
| < variable >_avg | HDF5 Dataset, 64-bit floating-point | optional |
| < variable >_rms | HDF5 Dataset, 64-bit floating-point | optional |
| < variable >_rey | HDF5 Dataset, 64-bit floating-point | optional |

Attributes:
|   |   |   |
|---|---|---|
| < variable >_wgt | 64-bit floating-point | optional/required |

If a variable exists, then the attribute is required.

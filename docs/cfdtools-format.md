# CFDtools internal formats

CFDtools defines several internal classes that aim at manipulating various types grids and data. One can currently find

- `meshbase.Mesh`
- `data.dataSet`
- `data.dataSetList`
- `vtkMesh`
- `vtkMeshList`

## mesh

A mesh is described by the positions of *nodes* and a connectivity between nodes. This connectivity can be based on cells or faces. In addition, at least boundaries should be tagged with a _physical names_.

### cell connectivity

### face connectivity

### tags

## data

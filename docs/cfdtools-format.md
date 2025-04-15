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

- `nfa_b` is the number of actual boundary faces
- `nfa_bp` is the number of periodic boundary faces
- `fa_count` is the total number of faces

### marks

Marks or tags are specified with `submeshmark` objects. Such mark contains the following properties:

- `name`
- `geodim`: node, face, or cell, possibly internal or bounding
- `type`: internal or bounding, can also define periodicity
- `index`: indirection table of elements (`indexlist` object)
- `connection` is a `meshconnection` object, need for mesh connections or periodic connections

### mesh connections


## data

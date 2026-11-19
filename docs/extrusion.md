# Extrusion

Extrusion transforms a 1D/2D mesh into a 2D/3D mesh by replicating it along a given direction,
stacking planes of the 1D/2D section to build a surface/volume mesh.

## Supported element types

| source element | target element | Order |
|------------|------------|-------|
| quad4      | hexa8      | linear |
| tri3       | penta6     | linear |
| quad9      | hexa27     | quadratic |

For quadratic elements, a mid-plane is inserted between each pair of planes to host
the edge-midpoint and face-center nodes of the 3D element.

## Python API

```python
import numpy as np

import cfdtools.cgns as cgns

input = cgns.cgnsMesh("file.cgns")
input.read_data()
mesh2d = input.export_mesh()

# Extrude 10 cells along z, uniformly spaced
mesh3d = mesh2d.export_extruded(
    direction=[0.0, 0.0, 1.0],
    extrude=np.linspace(0.0, 1.0, 11),  # 11 planes -> 10 cells
)
```

### Parameters

- **direction** (list, default `[0, 0, 1]`): unit vector along which to extrude.
- **extrude** (list of float, default `[0.0, 1.0]`): positions of each plane along
  `direction`. The number of planes determines the number of cells (`nplanes - 1`).
  For quadratic meshes, mid-plane nodes are placed at the midpoint of consecutive
  planes.
- **domain** (str, default `"fluid"`): prefix for the boundary marks created on the
  first and last planes (`domain + "0"` and `domain + "1"`).

### Boundary marks

Existing node-based boundary marks are replicated on every plane. Two additional
marks are created:

- `{domain}0` — nodes of the first (inlet) plane
- `{domain}1` — nodes of the last (outlet) plane

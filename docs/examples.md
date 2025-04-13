# Examples

## Readers

### From Legacy v3 IC3 Format to v4 IC3 (HDF5) Format

```python
import cfdtools.ic3.writerV4 as ic3writer
import cfdtools.ic3.reader_legacy as ic3reader

reader = ic3reader.reader("restart.000500.out")
reader.read_data()
rmesh = reader.export_mesh()

ic3write = ic3writer.writer(rmesh)
ic3write.write_data("mesh.h5", "solution.000000.h5")
```

## Writers

## Diagnosis

## Transformation

- [morphing/deformation](morph-sphere): internal generation of cube,morphed to a sphere and displayed as vtk
- [z-averaged or fourier modes](zconvolution) an extruded mesh

## Data

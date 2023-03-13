try:
    import pyvista as pv
    from pyvista import CellType
    importpyvista=True
except:
    importpyvista=False
import cfdtools.meshbase._elements as ele


map_ele = {
    'hexa8' : CellType.HEXAHEDRON
}

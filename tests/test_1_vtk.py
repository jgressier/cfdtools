import cfdtools.meshbase.simple as sm
from cfdtools.vtk import vtkMesh, vtkList
from cfdtools.hdf5 import h5File

import numpy as np
from pathlib import Path
import pytest


def test_cube_vtk(builddir):
    filename = "cube.vtu"
    filepath = builddir / filename
    cube = sm.Cube(10, 10, 10)
    mesh = cube.export_mesh()
    vtkmesh = vtkMesh(mesh)
    vtkmesh.write_data(filepath)
    Path(filepath).unlink()


def test_vtkread(datadir):
    name = datadir / "cubemixed0000.vtu"
    vtkfile = vtkMesh()
    vtkfile.read(name)
    assert vtkfile.pyvista_grid.n_cells == 1000
    assert vtkfile.volumes().sum() == pytest.approx(1.0)


def test_vtkdump(datadir, builddir):
    vtkname = datadir / "cubemixed0000.vtu"
    hname = builddir / "cubemixed0000.cfdh5"
    vtkfile = vtkMesh()
    vtkfile.read(vtkname)
    name = vtkfile.dumphdf(hname, overwrite=True)
    assert name == str(hname)  # since overwrite
    h5file = h5File(name)
    h5file.open()
    # assert h5file.datatype in ('unsvtk', 'dataset')
    assert h5file.datatype in ('dataset')
    assert "mesh" in h5file["/"].keys()
    vtkfile.importhdfgroup(h5file["mesh"])
    assert vtkfile.pyvista_grid.n_cells == 1000
    assert vtkfile.volumes().sum() == pytest.approx(1.0)

class Test_vtk_zconvolution():
    """Test the vtk_zconvolution function."""
    @pytest.fixture(autouse=True)
    def setup(self):
        cube = sm.Cube(10, 10, 10)
        mesh = cube.export_mesh()
        self.vtkmesh = vtkMesh(mesh)
        def var_u(xyz):
            dp = 2*np.pi
            return 10+2*np.sin(dp*xyz[:,2])+np.sin(2*dp*xyz[:,0])
        def var_p(xyz):
            dp = 2*np.pi
            return 100+np.sin(4*dp*xyz[:,2])*np.sin(2*dp*xyz[:,2])
        xyz = self.vtkmesh.pyvista_grid.cell_centers().points
        self.vtkmesh.pyvista_grid.cell_data['U'] = var_u(xyz)
        self.vtkmesh.pyvista_grid.cell_data['P'] = var_p(xyz)

    def test_default(self):
        assert self.vtkmesh.ncell == 1000
        vtkslice = self.vtkmesh.vtk_zconvolution()
        assert vtkslice.ncell == 100
        assert len(vtkslice.pyvista_grid.point_data.keys()) == 0
        assert len(vtkslice.pyvista_grid.cell_data.keys()) == 2
        for name in ('U_avg', 'P_avg'):
            assert name in vtkslice.pyvista_grid.cell_data.keys()

    def test_rms(self):
        assert self.vtkmesh.ncell == 1000
        vtkslice = self.vtkmesh.vtk_zconvolution(rms=True)
        assert vtkslice.ncell == 100
        assert len(vtkslice.pyvista_grid.point_data.keys()) == 0
        assert len(vtkslice.pyvista_grid.cell_data.keys()) == 4
        for name in ('U_avg', 'U_rms', 'P_avg', 'P_rms'):
            assert name in vtkslice.pyvista_grid.cell_data.keys()

    def test_snapshot(self):
        assert self.vtkmesh.ncell == 1000
        vtkslice = self.vtkmesh.vtk_zconvolution(snapshot=True)
        assert vtkslice.ncell == 100
        assert len(vtkslice.pyvista_grid.point_data.keys()) == 0
        assert len(vtkslice.pyvista_grid.cell_data.keys()) == 4
        for name in ('U_avg', 'U', 'P_avg', 'P'):
            assert name in vtkslice.pyvista_grid.cell_data.keys()

    def test_fourier(self):
        assert self.vtkmesh.ncell == 1000
        vtkslice = self.vtkmesh.vtk_zconvolution(nmode=2, phase=True)
        assert vtkslice.ncell == 100
        assert len(vtkslice.pyvista_grid.point_data.keys()) == 0
        assert len(vtkslice.pyvista_grid.cell_data.keys()) == 10
        for name in ('U_avg', 'P_avg'):
            assert name in vtkslice.pyvista_grid.cell_data.keys()
        vtkslice.write_data("test_fourier.vtu")

    def test_select(self):
        assert self.vtkmesh.ncell == 1000
        vtkslice = self.vtkmesh.vtk_zconvolution(nmode=1, rms=True, snapshot=True,
                                              select={'U': 'no_avg', 'P': [{'nmode': 2}, 'phase'] })
        assert vtkslice.ncell == 100
        assert len(vtkslice.pyvista_grid.point_data.keys()) == 0
        assert len(vtkslice.pyvista_grid.cell_data.keys()) == 10
        for name in ('U', 'U_k1', 'U_rms', 'P', 'P_avg', 'P_rms', 'P_k1', 'P_k2', 'P_k1_phase', 'P_k2_phase'):
            assert name in vtkslice.pyvista_grid.cell_data.keys()

def test_vtkList(datadir):
    namelist = list(datadir.glob("cubemixed00*.vtu"))
    vtklist = vtkList(namelist)
    assert vtklist.nfile == len(namelist)
    assert vtklist.allexist()


def test_vtkList_check(datadir):
    namelist = list(datadir.glob("cubemixed00*.vtu"))
    vtklist = vtkList(namelist, verbose=True)
    assert not vtklist.check_order('cellcenter')  # known to be mixed


def test_vtkList_read(datadir):
    namelist = list(datadir.glob("cubemixed00*.vtu"))
    vtklist = vtkList(namelist, verbose=True)
    vtklist.read()


def test_vtkList_dump(datadir, builddir):
    namelist = list(datadir.glob("cubemixed00*.vtu"))
    vtklist = vtkList(namelist, verbose=True)
    vtklist.read()
    h5filename = builddir / "vtklist.hdf"
    newname = vtklist.dumphdf(h5filename)
    h5file = h5File(newname)
    h5file.open()
    assert h5file.datatype == 'datalist'
    assert len(h5file["/datalist"].keys()) == 10
    Path(newname).unlink()


def test_vtkList_dumpxdmf(datadir, tmpdir):
    """Test the xmf file creation."""
    namelist = sorted(list(datadir.glob("cubemixed000[0-1].vtu")))
    vtklist = vtkList(namelist, verbose=True)
    vtklist.read()
    # vtu files do not contain time values, add them artificially
    for dataset in vtklist._data._datalist:
        dataset['time'] = np.array([3.14])

    h5filename = tmpdir / "vtklist.hdf"
    vtklist.dumphdf(h5filename, xdmf=True)

    xmf_filepath = tmpdir / "vtklist.xmf"
    assert xmf_filepath.exists()

    xdmf_content_formatted = """
        <Xdmf Version="3.0">
            <Domain>
                <Grid Name="IC3" GridType="Collection" CollectionType="Temporal">
        """
    for i in range(2):
        xdmf_content_formatted += f"""
                    <Grid Name="Unstructured Mesh">
                        <Time Value="{i}"/>
                        <Geometry GeometryType="XYZ">
                            <DataItem Dimensions="3993" Format="HDF">vtklist.hdf:/mesh/nodes</DataItem>
                        </Geometry>
                        <Topology NumberOfElements="1000" TopologyType="Hexahedron">
                            <DataItem Dimensions="8000" Format="HDF">vtklist.hdf:/mesh/cells/hexa8</DataItem>
                        </Topology>
                        <Attribute AttributeType="Scalar" Center="Cell" Name="Q">
                            <DataItem Dimensions="1000" Format="HDF">vtklist.hdf:/datalist/i{i:06}/Q</DataItem>
                        </Attribute>
                    </Grid>
        """ 
    xdmf_content_formatted += """
                </Grid>
            </Domain>
        </Xdmf>
        """
    xdmf_content = ''.join([t.strip() for t in xdmf_content_formatted.split('\n')])

    lines = open(xmf_filepath, 'r').read()
    # remove end of lines
    lines = lines.replace('\n', '')
    # replace temporary path with reference path
    lines = lines.replace(str(h5filename), "vtklist.hdf")
    assert xdmf_content == lines

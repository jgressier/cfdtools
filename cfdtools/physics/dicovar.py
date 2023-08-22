# -*- coding: utf-8 -*-
"""
Description
-----------

Module containing the definition of a DicoVar object, which allow the switch
from different variable denomination in vtk output files, depending on which
solver denomination you use (CharlesX, OpenFOAM, Hybrid, ...)

Contains
--------
dicoVar : _DicoVar object
    The current database for variables denomination
setDicoVar(software) : function
    Set the current database to a given software denomination
printDatabase() : function
    Print the current database denomination

Examples
--------
When implementing a vtk file manipulation function

>>> # Import dicovar module database
>>> from hades.common.dicovar import dicoVar as dv
>>> # test if the field U_AVG is present
>>> if data_outVTK.GetPointData().HasArray(dv('u_avg')) != 1:
>>>     raise ValueError("Error : field U_AVG not present")
>>> # test if the field RHO_AVG is present
>>> if data_outVTK.GetPointData().HasArray(dv['rho_avg']) != 1:
>>>     raise ValueError("Error : field RHO_AVG not present")
>>> # test if the field TauWallAvg is present
>>> if wall_outVTK.GetPointData().HasArray(dv.get('tauw_avg')) != 1:
        raise ValueError("Error : field TauWallAvg not present")

Notes
-----
- The implementation allow to use several names to set one software (ex : 'cx',
  'CharlesX', ... for CharlesX denomination)
- Several ways to get a variables : dv[...], dv(...), dv.get(...)
- variable name handles capital letters : dv['U_avg'] <=> dv['u_avg']
"""
import copy as _copy

# -----------------------------------------------------------------
# Part to modify in order to implement other variable denominations
# ---------------

# CharlesX variables dico
_cxNames = ['IC3', 'ic3', 'CharlesX', 'cx', 'charlesx']
_cxDico = {
    'P': 'P',
    'p': 'P',
    'Ps': 'P',
    'rho': 'RHO',
    'T': 'T',
    'V': 'U',
    'Vx': 'U_X',
    'Vy': 'U_Y',
    'Vz': 'U_Z',
    'u_avg': 'U_AVG',
    'rho_avg': 'RHO_AVG',
    'tauw_avg': 'TauWallAvg',
    'mulam_avg': 'MU_LAM_AVG',
    't_avg': 'T_AVG',
    'u_rms': 'U_RMS',
    'u_rey': 'U_REY',
}

_default = 'IC3'

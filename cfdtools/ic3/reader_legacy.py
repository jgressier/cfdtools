# coding: utf8

# Import modules
import struct
import numpy as np
import sys
import cfdtools.api as api
import cfdtools.meshbase._mesh as _mesh
import cfdtools.meshbase._data as _data
#import cfdtools.ic3._ic3 as ic3b
from cfdtools.ic3._ic3 import type2nbytes, restartSectionHeader, ic3_restart_codes, BinaryRead, zonekind2type, nno2fatype, properties_ugpcode

###################################################################################################

@api.fileformat_reader('IC3', '.ic3')
class reader(api._files):
    '''Implementation of the reader to read IC3 restart files.'''

    def __init__(self, filename, cIntegrity=False):
        '''
        Initialization of an IC3 restart reader.
        Just save the filename and a boolean for an integrity check.
        input   : IC3 restart file name [type string]
                  whether to check the integrity of the file beforehand [type boolean]
        '''
        super().__init__(filename)
        self.check_integrity = cIntegrity
        self.ic3_version = -1

    def __str__(self):
        s = '    filename: '+self.filename
        s+= '\n   simulation: '+str(list(self.simulation_state.keys()))
        s+= '\n    mesh keys: '+str(list(self.mesh.keys()))
        s+= '\nvariable keys: '+str(list(self.variables.keys()))
        return s

    def printinfo(self):
        print(self)
        print("- mesh properties")
        for key,item in self.mesh.items():
            if isinstance(item,dict):
                for key2,item2 in item.items():
                    if isinstance(item2,dict):
                        for key3,item3 in item2.items():
                            api._printreadable('  mesh.'+key+'.'+key2+'.'+key3, item3)
                    else:
                        api._printreadable('  mesh.'+key+'.'+key2, item2)
            else:
                api._printreadable('  mesh.'+key, item)
        print("- variable properties")
        for key,item in self.variables.items():
            for key2,value in item.items():
                api._printreadable('variables.'+key+'.'+key2, value)
                # for key3,value3 in value.items():
                #     print('variables.'+key+'.'+key2+'.'+key3+':'+str(type(value3)))

    def read_data(self):
        '''
        Main method of the IC3 restart reader.
        Parses in order the file using sub-methods described below.
        output  : the mesh itself
                  the lot of variables stored in the restart file
                  information on the state of the simulation
        '''
        api.io.print('std',":: READER RESTART IC3 ::")

        if not self._exists:
            print("Fatal error. File %s cannot be found."%(self.filename))
            exit()

        # Open the file for binary reading
        api.io.print('debug','reading ',self.filename)
        self.fid = open(self.filename, "rb")

        api.io.print('std', "Reading header ..")
        self.__ReadRestartHeader()
        #
        api.io.print('std', "Reading connectivity ..")
        self.__ReadRestartConnectivity()
        #
        api.io.print('std', "Reading informative values ..")
        self.__ReadInformativeValues()
        #
        api.io.print('std', "Reading variables ..")
        self.__ReadRestartVar()
        #

        # Before returning, close the file
        self.fid.close()
        del self.fid

        #return self.mesh["coordinates"], self.mesh["connectivity"]["e2v"], self.mesh["bocos"], self.variables["nodes"], self.variables["cells"], (self.simulation_state, self.mesh["params"])
        meshdata = _mesh.mesh(self.mesh['params']['cv_count'], self.mesh['params']['no_count'])
        meshdata.set_nodescoord_nd(self.mesh['coordinates'])
        meshdata.set_face2cell(self.mesh['connectivity']['cvofa'])
        meshdata.set_face2node(self.mesh['connectivity']['noofa'])
        meshdata.set_bocos(self.mesh['bocos'])
        meshdata.set_celldata(self.variables['cells'])
        meshdata.set_nodedata(self.variables['nodes'])
        meshdata.set_facedata(self.variables['faces'])
        meshdata.set_params(self.mesh['params'])
        meshdata.update_params(self.simulation_state)
        if 'partition' in self.mesh.keys():
            meshdata.set_partition(self.mesh['partition'])

        return meshdata

    def __ReadRestartHeader(self):
        '''
        Method reading the header of a restart file.
        It is composed of two integers, the "magic number" used as a flag for endianness
        and the CharlesX version number.
        input   : handle on an open restart file, [type file identifier]
        output  : the endianness of the open restart file [type boolean]
        '''

        # By default suppose big-endian format
        self.byte_swap = False

        # Read the first integer (int64)
        s = list(BinaryRead(self.fid, "ii", False, 2*type2nbytes["int32"]))
        # If, with big-endian assumption, the first integer comes out wrong, swap to little-endian
        if s[0] != ic3_restart_codes["UGP_IO_MAGIC_NUMBER"]:
            # Change the flag
            self.byte_swap=True
            # Transform the second integer of the list to match the version number
            aux_struct = struct.Struct(">i")
            packed_version = aux_struct.pack(s[1])
            del aux_struct
            aux_struct = struct.Struct("<i")
            s[1] = aux_struct.unpack(packed_version)[0]

        # Some info for the user
        self.ic3_version = s[1]
        api.io.print('std', f"\t Restart file is version {self.ic3_version}")

    def __ReadRestartConnectivity(self):
        '''
        Method reading the first blocks passed the header, containing informations
        on the nodes, the faces, the cells and the connectivity in between those
        input   : handle on an open restart file [type file identifier]
                  endianness flag [boolean]
        output  : mesh structure containing all the necessary information to build the grid
        '''

        # Initialize the mesh
        self.mesh = {"params":{"no_count":0, "fa_count":0, "cv_count":0, "noofa_count":0, "nboco":0},
                     "connectivity":{"noofa":{}, "cvofa":{}, "nkeys":0},
                     "coordinates":None,
                     "bocos":{"nfa_b":0, "nfa_bp":0},
                     "partition":None}

        # Get the header
        h = restartSectionHeader()
        if (not h.readVar(self.fid, self.byte_swap,["UGP_IO_NO_FA_CV_NOOFA_COUNTS"])): exit()

        # Store the size informations at the right place
        self.mesh["params"]["no_count"] = h.idata[0]
        self.mesh["params"]["fa_count"] = h.idata[1]
        self.mesh["params"]["cv_count"] = h.idata[2]
        self.mesh["params"]["noofa_count"] = h.idata[3]
        api.io.print("std", "mesh with {} cells, {} faces and {} ndoes".format(h.idata[2], h.idata[1], h.idata[0]))
        del h

        # Integrity check
        if self.check_integrity:
            # Check the restart is whole by counting the global ids of no, fa and cv
            # For nodes
            api.io.print('std', "\t Checking nodes integrity .."); sys.stdout.flush()
            h = restartSectionHeader()
            if(not h.readVar(self.fid, self.byte_swap,["UGP_IO_NO_CHECK"])): exit()

            assert h.idata[0] == self.mesh["params"]["no_count"]
            assert h.id[0] == ic3_restart_codes["UGP_IO_NO_CHECK"]
            nodes_id = np.empty((self.mesh["params"]["no_count"],), dtype=np.int32)
            for loopi in range(self.mesh["params"]["no_count"]):  # xrange to range (python3 portage)
                s = BinaryRead(self.fid, "i", self.byte_swap, type2nbytes["int32"])
                nodes_id[loopi] = s[0]
            assert np.allclose(nodes_id, np.arange(self.mesh["params"]["no_count"]))
            api.io.print('std', "ok.") ; sys.stdout.flush()
            del nodes_id, h
            # For faces
            api.io.print('std', "\t Checking faces integrity .."); sys.stdout.flush()
            h = restartSectionHeader()
            if(not h.readVar(self.fid, self.byte_swap,["UGP_IO_FA_CHECK"])): exit()
            faces_id = np.empty((self.mesh["params"]["fa_count"],), dtype=np.int32)
            for loopi in range(self.mesh["params"]["fa_count"]):  # xrange to range (python3 portage)
                s = BinaryRead(self.fid, "i", self.byte_swap, type2nbytes["int32"])
                faces_id[loopi] = s[0]
            assert np.allclose(faces_id, np.arange(self.mesh["params"]["fa_count"]))
            api.io.print('std', "ok."); sys.stdout.flush()
            del faces_id, h
            # For cells
            api.io.print('std', "\t Checking cells integrity .."); sys.stdout.flush()
            h = restartSectionHeader()
            if(not h.readVar(self.fid, self.byte_swap,["UGP_IO_CV_CHECK"])): exit()
            assert h.idata[0] == self.mesh["params"]["cv_count"]
            assert h.id[0] == ic3_restart_codes["UGP_IO_CV_CHECK"]
            cells_id = np.empty((self.mesh["params"]["cv_count"],), dtype=np.int32)
            for loopi in range(self.mesh["params"]["cv_count"]):  # xrange to range (python3 portage)
                s = BinaryRead(self.fid, "i", self.byte_swap, type2nbytes["int32"])
                cells_id[loopi] = s[0]
            assert np.allclose(cells_id, np.arange(self.mesh["params"]["cv_count"]))
            api.io.print('std', "ok."); sys.stdout.flush()
            del cells_id, h

        # The two connectivities now
        #
        #- First, NOOFA
        #
        api.io.print('std', "\t Parsing face to node connectivity .."); sys.stdout.flush()
        h = restartSectionHeader()
        if (not h.readVar(self.fid, self.byte_swap,["UGP_IO_NOOFA_I_AND_V"])): exit()
        #
        assert h.idata[0] == self.mesh["params"]["fa_count"]
        assert h.idata[1] == self.mesh["params"]["noofa_count"]
        # Get the node count per face
        nno_per_face = np.empty((self.mesh["params"]["fa_count"],), dtype=np.int32)
        for loopi in range(self.mesh["params"]["fa_count"]): # xrange to range (python3 portage)
            s = BinaryRead(self.fid, "i", self.byte_swap, type2nbytes["int32"])
            nno_per_face[loopi] = s[0]
        uniq, counts = np.unique(nno_per_face, return_counts=True)
        #print("uniq:",uniq)
        uniq = [nno2fatype[val] for val in uniq]
        api.io.print('std', "found %s faces .."%(" and ".join(uniq))); sys.stdout.flush()
        # Initialize the proper connectivity arrays in self.mesh
        self.mesh["connectivity"]["noofa"]["listofStarts_f2v"] = np.concatenate(([0,], np.cumsum(nno_per_face)))
        self.mesh["connectivity"]["noofa"]["face2vertex"] = np.zeros((np.sum(nno_per_face),), dtype=np.int64)
        assert self.mesh["connectivity"]["noofa"]["listofStarts_f2v"][0] == 0
        assert self.mesh["connectivity"]["noofa"]["listofStarts_f2v"][-2] == self.mesh["connectivity"]["noofa"]["face2vertex"].size - nno_per_face[-1]
        # Now loop on the restart file to fill the connectivities
        for loopi in range(self.mesh["params"]["fa_count"]):  # xrange to range (python3 portage)
            sta, sto = self.mesh["connectivity"]["noofa"]["listofStarts_f2v"][loopi], self.mesh["connectivity"]["noofa"]["listofStarts_f2v"][loopi+1]
            s = BinaryRead(self.fid, "i"*nno_per_face[loopi], self.byte_swap, type2nbytes["int32"]*nno_per_face[loopi])
            self.mesh["connectivity"]["noofa"]["face2vertex"][sta:sto] = np.asarray(s).astype(np.int64)
        api.io.print('std', "ok."); sys.stdout.flush()
        del nno_per_face, h, uniq, counts
        #
        #- Second, CVOFA
        #
        api.io.print('std', "\t Parsing face to cell connectivity .."); sys.stdout.flush()
        h = restartSectionHeader()
        if (not h.readVar(self.fid, self.byte_swap,["UGP_IO_CVOFA"])): exit()

        assert h.idata[0] == self.mesh["params"]["fa_count"]
        assert h.idata[1] == 2
        # Initialize the proper connectivity arrays in self.mesh
        self.mesh["connectivity"]["cvofa"]["cvofa"] = np.zeros((self.mesh["params"]["fa_count"], 2), dtype=np.int64)
        # Now loop on the restart file to fill the connectivities
        for loopi in range(self.mesh["params"]["fa_count"]):  # xrange to range (python3 portage)
            s = BinaryRead(self.fid, "ii", self.byte_swap, type2nbytes["int32"]*self.mesh["connectivity"]["cvofa"]["cvofa"].shape[1])
            self.mesh["connectivity"]["cvofa"]["cvofa"][loopi, :] = np.asarray(s).astype(np.int64)
        api.io.print('std', "ok."); sys.stdout.flush()
        #print("RA",self.mesh["connectivity"]["cvofa"]["cvofa"])
        del h
        # Checks and a few associations
        assert self.mesh["connectivity"]["cvofa"]["cvofa"].max() < self.mesh["params"]["cv_count"]
        assert self.mesh["connectivity"]["cvofa"]["cvofa"].max() == self.mesh["params"]["cv_count"]-1
        uniq, counts = np.unique(self.mesh["connectivity"]["cvofa"]["cvofa"][:,1], return_counts=True)
        # Number of assigned boundary faces
        try:
            iwhere = np.where(uniq==-1)[0][0]
        except IndexError:
            self.mesh["bocos"]["nfa_b"] = 0
        else:
            self.mesh["bocos"]["nfa_b"] = counts[iwhere]
        del uniq, counts
        # Number of periodic boundary faces
        self.mesh["bocos"]["nfa_bp"] = np.count_nonzero(self.mesh["connectivity"]["cvofa"]["cvofa"][:,1] < -1)
        # All periodic boundary faces to -1
        # mask = self.mesh["connectivity"]["cvofa"]["cvofa"][:,1] < -1
        # self.mesh["connectivity"]["cvofa"]["cvofa"][:,1][mask] = -1
        #print("RC",self.mesh["connectivity"]["cvofa"]["cvofa"])
        #
        # The boundary conditions now
        api.io.print('std', "\t Parsing boundary conditions .."); sys.stdout.flush()
        while True:
            h = restartSectionHeader()
            if(not h.readVar(self.fid, self.byte_swap,["UGP_IO_FA_ZONE"],reset_offset=False)): break

            self.mesh["params"]["nboco"] += 1
            self.mesh["bocos"][h.name] = {}
            self.mesh["bocos"][h.name]["type"] = zonekind2type[h.idata[0]]
            self.mesh["bocos"][h.name]["fa_range"] = np.array([h.idata[1], h.idata[2]])
            self.mesh["bocos"][h.name]["periodic_transform"] = h.rdata
            #
            famin, famax = self.mesh["bocos"][h.name]["fa_range"]
            sta = self.mesh["connectivity"]["noofa"]["listofStarts_f2v"][famin]
            try:
                sto = self.mesh["connectivity"]["noofa"]["listofStarts_f2v"][famax+1]
            except IndexError:
                sto = self.mesh["connectivity"]["noofa"]["face2vertex"].size
            #
            self.mesh["bocos"][h.name]["slicing"] = np.unique(self.mesh["connectivity"]["noofa"]["face2vertex"][sta:sto])
            if h.idata[0] == 6:
                break
        api.io.print('standard', "ok.")
        sys.stdout.flush()

        # Parse the header of the partition information
        api.io.print('std', "\t Parsing partitioning information ..."); sys.stdout.flush()
        h = restartSectionHeader()
        if(not h.readVar(self.fid, self.byte_swap,["UGP_IO_CV_PART"])): exit()

        self.mesh["partition"] = {}
        self.mesh["partition"]['npart'] = h.idata[1]
        self.mesh["partition"]['icvpart'] = np.zeros((self.mesh["params"]["cv_count"],), dtype=np.int32)
        for loopi in range(self.mesh["params"]["cv_count"]):  # xrange to range (python3 portage)
            s = BinaryRead(self.fid, "i", self.byte_swap, type2nbytes["int32"])
            self.mesh["partition"]['icvpart'][loopi] = s[0]
        #print(h)
        api.io.print('std', "ok.")
        sys.stdout.flush()

        # The coordinates of the vertices finally
        api.io.print('std', "\t Parsing vertices coordinates .."); sys.stdout.flush()
        h = restartSectionHeader()
        if(not h.readVar(self.fid, self.byte_swap,["UGP_IO_X_NO"])): exit()

        assert h.idata[0] == self.mesh["params"]["no_count"]
        assert h.idata[1] == 3
        self.mesh["coordinates"] = np.zeros((self.mesh["params"]["no_count"], 3), dtype=np.float64)
        for loopi in range(self.mesh["params"]["no_count"]):  # xrange to range (python3 portage)
            s = BinaryRead(self.fid, "ddd", self.byte_swap, type2nbytes["float64"]*self.mesh["coordinates"].shape[1])
            self.mesh["coordinates"][loopi, :] = np.asarray(s)
        self.mesh["coordinates"] = np.ascontiguousarray(self.mesh["coordinates"])
        api.io.print('std', "ok."); sys.stdout.flush()

    def __ReadInformativeValues(self):
        '''
        Method reading all the values also stored in a restart file,
        i.e. the step number, the time, the timestep.
        input   : handle on an open restart file, [type file identifier]
                  endianness flag [boolean]
        output  : simulation state structure containing informations about the current state
                  of the simulation
        '''

        # Initialize the state dictionary
        self.simulation_state = {"step":0, "dt":0, "time":0, "wgt":{}}

        # First, a header saying data to introduce to this block we are parsing now


        #removed as not used at the moment
        #h = restartSectionHeader()
        #if(not readVar(self.fid, self.byte_swap,"UGP_IO_DAT)A")
        #if (not varfound): exit()

        reset_offset=True
        while True:
            h = restartSectionHeader()
            if (not h.readVar(self.fid, self.byte_swap,["UGP_IO_I0"],reset_offset=reset_offset)): break
            else: reset_offset=False
            print("I0 var "+h.name.lower())
            self.simulation_state[h.name.lower()] = h.idata[0]
        reset_offset=True
        while True:
            h = restartSectionHeader()
            if (not h.readVar(self.fid, self.byte_swap,["UGP_IO_D0"],reset_offset=reset_offset)): break
            else: reset_offset=False
            print("D0 var "+h.name.lower())
            self.simulation_state[h.name.lower()] = h.rdata[0]

    # def get_datacell_properties(self):
    #     return self.variables["cells"]["_info"]

    # def get_ndof(self):
    #     return self.get_datacell_properties["ndof"]

    def _set_ndof_properties(self, intndof):
        ndof = 1
        if self.ic3_version < 0:
            raise ValueError("unknown IC3 version number")
        elif self.ic3_version < 3:
            pass
            # ignore ints 
            # if intndof != 0: # 0 is expected value for version 1 and 2
            #     api.io.print("warning", "unexpected non zero value for ndof in IC3 version 1 and 2")
        else: # version >= 3
            if intndof == 0: # 0 is NOT expected 
                api.io.print("warning", "unexpected zero value for ndof in IC3 version 3; set to 1 !")
                ndof = 1
            else:
                ndof = intndof
        # if self.get_datacell_properties().get("ndof", ndof) != ndof:
        #     raise ValueError("Inconsistent cell data size (ndof)")
        # self.get_datacell_properties()["ndof"] = ndof
        return ndof

    def __ReadRestartVar(self):
        '''
        Method reading the variables from the restart file
        input   : handle on an open restart file, [type file identifier]
                  endianness flag [boolean]
                  mesh structure
        output  : structure containing all the variables
        '''

        # Some extra modules
        import copy

        # Initialize the variable dictionary
        ncv = self.mesh["params"]["cv_count"]
        self.variables = {"nodes":{}, "cells":{}, "faces":{}}

        # First come the scalars
        api.io.print('std', "\t First the scalars ...")
        reset_offset=True
        while True:
            h = restartSectionHeader()
            if (not h.readVar(self.fid, self.byte_swap,
                ["UGP_IO_NO_D1","UGP_IO_NO_II1","UGP_IO_FA_D1","UGP_IO_CV_D1","UGP_IO_CV_II1"],
                reset_offset=reset_offset)): break
            reset_offset=False
            #
            typechar = properties_ugpcode[h.id[0]]['structcode']
            typesize = properties_ugpcode[h.id[0]]['size']
            nptype = properties_ugpcode[h.id[0]]['numpytype']
            #
            if h.idata[0] == self.mesh["params"]["no_count"]:
                self.variables["nodes"][h.name] = np.zeros((self.mesh["params"]["no_count"],), dtype=nptype)
                s = BinaryRead(self.fid, typechar*self.mesh["params"]["no_count"], self.byte_swap, typesize*self.mesh["params"]["no_count"])
                self.variables["nodes"][h.name] = np.asarray(s)
                api.io.print('std', "\t %s%s:\t %+.5e / %+.5e / %+.5e (min/mean/max)."%(h.name, ' '*(20-len(h.name)), np.asarray(s).min(), np.mean(np.asarray(s)), np.asarray(s).max()))
            elif h.idata[0] == self.mesh["params"]["fa_count"]:
                self.variables["faces"][h.name] = np.zeros((self.mesh["params"]["fa_count"],), dtype=nptype)
                s = BinaryRead(self.fid, typechar*self.mesh["params"]["fa_count"], self.byte_swap, typesize*self.mesh["params"]["fa_count"])
                self.variables["faces"][h.name] = np.asarray(s)
                api.io.print('std', "\t %s%s:\t %+.5e / %+.5e / %+.5e (min/mean/max)."%(h.name, ' '*(20-len(h.name)), np.asarray(s).min(), np.mean(np.asarray(s)), np.asarray(s).max()))
            elif h.idata[0] == ncv:
                api.io.print("internal", "cell variable section of size {}x{}".format(h.idata[0], h.idata[1]))
                ndof = self._set_ndof_properties(h.idata[1])
                pdata = np.zeros((ndof*ncv,), dtype=nptype)
                s = BinaryRead(self.fid, typechar*ndof*ncv, self.byte_swap, typesize*ndof*ncv)
                pdata = np.asarray(s)
                # If multiple connectivities, gotta order the tables correctly
                if self.mesh["connectivity"]["nkeys"] > 1:
                    aux = copy.deepcopy(pdata)
                    pdata = np.empty((0,), dtype=aux.dtype)
                    for uns_type, indices in self.mesh["connectivity"]["cell_indices"].items():
                        pdata = np.concatenate((pdata, aux[indices]))
                self.variables["cells"][h.name] = _data.celldata(ndof=ndof)
                self.variables["cells"][h.name].set_data(pdata)
                api.io.print('std', "\t %s%s:\t %+.5e / %+.5e / %+.5e (min/mean/max)."%(h.name, ' '*(20-len(h.name)), np.asarray(s).min(), np.mean(np.asarray(s)), np.asarray(s).max()))
            else:
                api.io.print('std', "Fatal error. Incoherence in dataset %s. Exiting."%(h.name))
                exit()
        api.io.print('std', "\t ok.")
        

        # Then the vectors
        api.io.print('std', "\t Then the vectors ..")
        reset_offset=True
        while True:
            h = restartSectionHeader()
            if(not h.readVar(self.fid, self.byte_swap,["UGP_IO_NO_D3","UGP_IO_FA_D3","UGP_IO_CV_D3"],reset_offset=reset_offset)): break
            reset_offset=False

            if h.idata[0] == self.mesh["params"]["no_count"]:
                self.variables["nodes"][h.name] = np.zeros((self.mesh["params"]["no_count"], 3), dtype=np.float64)
                s = BinaryRead(self.fid, "ddd"*self.mesh["params"]["no_count"], self.byte_swap, type2nbytes["float64"]*self.mesh["params"]["no_count"]*3)
                self.variables["nodes"][h.name] = np.asarray(s).reshape((-1, 3))
                api.io.print('std', "\t %s%s:\t %+.5e / %+.5e / %+.5e (min/mean/max)."%(h.name, ' '*(20-len(h.name)), np.asarray(s).min(), np.mean(np.asarray(s)), np.asarray(s).max()))
            elif h.idata[0] == self.mesh["params"]["fa_count"]:
                self.variables["faces"][h.name] = np.zeros((self.mesh["params"]["fa_count"], 3), dtype=np.float64)
                s = BinaryRead(self.fid, "ddd"*self.mesh["params"]["fa_count"], self.byte_swap, type2nbytes["float64"]*self.mesh["params"]["fa_count"]*3)
                self.variables["faces"][h.name] = np.asarray(s).reshape((-1, 3))
                api.io.print('std', "\t %s%s:\t %+.5e / %+.5e / %+.5e (min/mean/max)."%(h.name, ' '*(20-len(h.name)), np.asarray(s).min(), np.mean(np.asarray(s)), np.asarray(s).max()))
            elif h.idata[0] == ncv:
                ndof = self._set_ndof_properties(h.idata[1])
                pdata = np.zeros((ndof*ncv, 3), dtype=np.float64)
                s = BinaryRead(self.fid, "ddd"*ndof*ncv, self.byte_swap, type2nbytes["float64"]*ndof*ncv*3)
                pdata = np.asarray(s).reshape((-1, 3))
                # If multiple connectivities, gotta order the tables correctly
                if self.mesh["connectivity"]["nkeys"] > 1:
                    aux = copy.deepcopy(pdata)
                    pdata = np.empty((0, 3), dtype=aux.dtype)
                    for uns_type, indices in self.mesh["connectivity"]["cell_indices"].items():
                        pdata = np.concatenate((pdata, aux[indices, :]), axis=0)
                self.variables["cells"][h.name] = _data.celldata(ndof=ndof)
                self.variables["cells"][h.name].set_data(pdata)
                api.io.print('std', "\t %s%s:\t %+.5e / %+.5e / %+.5e (min/mean/max)."%(h.name, ' '*(20-len(h.name)), np.asarray(s).min(), np.mean(np.asarray(s)), np.asarray(s).max()))
            else:
                api.io.print('std', "Fatal error. Incoherence in dataset %s. Exiting."%(h.name))
                exit()
        api.io.print('std', "\t ok.")
        

        # Then the tensors
        api.io.print('std', "\t Then the tensors ..")
        reset_offset=True
        while True:
            h = restartSectionHeader()
            if(not h.readVar(self.fid, self.byte_swap,["UGP_IO_CV_D33"],reset_offset=reset_offset)): break
            reset_offset=False

            if h.idata[0] == ncv:
                ndof = self._set_ndof_properties(h.idata[1])
                pdata = np.zeros((ndof*ncv, 3, 3), dtype=np.float64)
                s = BinaryRead(self.fid, "d"*ndof*ncv*3*3, self.byte_swap, type2nbytes["float64"]*ndof*ncv*3*3)
                pdata = np.asarray(s).reshape((-1, 3, 3))
                # If multiple connectivities, gotta order the tables correctly
                if self.mesh["connectivity"]["nkeys"] > 1:
                    aux = copy.deepcopy(pdata)
                    pdata = np.empty((0, 3, 3), dtype=aux.dtype)
                    for uns_type, indices in self.mesh["connectivity"]["cell_indices"].items():
                        pdata = np.concatenate((pdata, aux[indices, :, :]), axis=0)
                self.variables["cells"][h.name] = _data.celldata(ndof=ndof)
                self.variables["cells"][h.name].set_data(pdata)
                api.io.print('std', "\t %s%s:\t %+.5e / %+.5e / %+.5e (min/mean/max)."%(h.name, ' '*(20-len(h.name)), np.asarray(s).min(), np.mean(np.asarray(s)), np.asarray(s).max()))
            else:
                api.io.print('std', "Fatal error. Incoherence in dataset %s. Exiting."%(h.name))
                exit()
        api.io.print('std', "\t ok.")
        

    # def __reachedEOF(self):
    #     '''
    #     Method to check whether the reader reached EOF.
    #     It should happen right after reading all the variable blocks.
    #     input   : handle on an open restart file, [type file identifier]
    #               endianness flag [boolean]
    #     output  : check on end-of-file [boolean]
    #     '''

    #     # Read the next header
    #     h = restartSectionHeader()
    #     h.read(self.fid, self.byte_swap)
    #     if h.name != "EOF" or\
    #        h.id != ic3_restart_codes["UGP_IO_EOF"]:
    #        return False

    #     return True

###################################################################################################

# if __name__ == "__main__":
#     '''
#     The script is supposed to be used with command line arguments
#     but if it is not, it runs a test on a pre-defined file name.
#     '''

#     # Module import for script use
#     import copy
#     import os
#     import sys
#     import vtk as _vtk
#     from vtk.util import numpy_support

#     # Parse command line arguments seeking a restart file name
#     defaultFName = "restart_test.out"
#     if len(sys.argv) == 1:
#         fName = defaultFName
#         cIntegrity = False
#         api.io.print('std', "Use requires the name of a restart file to be input in argument.")
#         api.io.print('std', "Using a default name, %s."%(defaultFName))
#     else:
#         fName = defaultFName
#         cIntegrity = False
#         for arg in sys.argv[1:]:
#             if '--restartname' in arg:
#                 fName = arg.split('=')[-1]
#             elif '--checkIntegrity' in arg:
#                 _cIntegrity = arg.split('=')[-1]
#                 if _cIntegrity == "False":
#                     cIntegrity = False
#                 elif _cIntegrity == "True":
#                     cIntegrity = True
#     if fName == defaultFName:
#         api.io.print('std', "Use requires the name of a restart file to be input in argument.")
#         api.io.print('std', "Using a default name, %s."%(defaultFName))

#     # Check given file name
#     if not os.path.isfile(fName):
#         api.io.print('std', "Fatal error: %s cannot be found. Exiting."%(fName))
#         exit()

#     # Else proceed
#     readr = ReaderRestartIC3(fName, cIntegrity)
#     xyz, co, bocos, simulation_state, nodesvar, cellsvar, params = readr.read_data()

#     ###################################################################################################

#     # Now, we can export the data under a more human-friendly format
#     #
#     #- IN VTK
#     #
#     # Build point coordinates for VTK
#     nb_point = xyz.shape[0]
#     if hasattr(_vtk, 'vtkSOADataArrayTemplate'):
#         # VTK 8.1.0+: use SOA coordinate array data structure (zero-copy)
#         points_vtk = _vtk.vtkSOADataArrayTemplate[np.float64]()  # shall avoid copy
#         points_vtk.SetNumberOfComponents(3)
#         points_vtk.charlesx_arrays = []
#         for index in range(3):
#             np_array = copy.deepcopy(xyz[:, index])
#             points_vtk.SetArray(index, np_array, nb_point, True, True)
#             points_vtk.charlesx_arrays.append(np_array)  # record array ref so it won't be gc'ed
#     else:
#         # VTK 8.0.1-: build a new AOS coordinate array
#         points_vtk = numpy_support.numpy_to_vtk(mesh["coordinates"], deep=True)
#     vtkpoints = _vtk.vtkPoints()
#     vtkpoints.SetData(points_vtk)
#     # Build connectivity for VTK
#     vtkprimitive = {
#         'bi': _vtk.vtkLine(),
#         'tri': _vtk.vtkTriangle(),
#         'qua': _vtk.vtkQuad(),
#         'tet': _vtk.vtkTetra(),
#         'hex': _vtk.vtkHexahedron(),
#         'pri': _vtk.vtkWedge(),
#         'pyr': _vtk.vtkPyramid(),
#     }
#     cells = np.empty((0, ), dtype=np.int64)
#     cell_types = np.empty((0, ), dtype=np.int64)
#     offsets = np.empty((0, ), dtype=np.int64)
#     offset_start = 0
#     total_nb_cells = 0

#     for uns_type, connect in co.items():

#         nb_cells = connect.shape[0]

#         # record cell-type for each individual cell
#         if uns_type in vtkprimitive:
#             cell_type = vtkprimitive[uns_type].GetCellType()
#         else:
#             raise ValueError(str_error('Unknown cell type'))
#         cell_types = np.append(cell_types, np.tile(cell_type, (nb_cells, 1)))

#         # put number of vertices before each cell
#         cells = np.append(cells, np.concatenate((np.tile(nodes_per_cell[uns_type], (nb_cells, 1)), connect),
#                                                 axis=1).flat)

#         # start offset of each cell in connectivity array
#         offset_stop = offset_start + nb_cells * (nodes_per_cell[uns_type] + 1)
#         offsets = np.append(offsets, np.arange(offset_start, offset_stop, nodes_per_cell[uns_type] + 1))
#         offset_start = offset_stop
#         total_nb_cells += nb_cells
#     idtype_vtk = _vtk.vtkIdTypeArray().GetDataType()
#     uchartype_vtk = _vtk.vtkUnsignedCharArray().GetDataType()

#     cells_vtk = numpy_support.numpy_to_vtk(cells, deep=True, array_type=idtype_vtk)
#     cell_array = _vtk.vtkCellArray()
#     cell_array.SetCells(total_nb_cells, cells_vtk)
#     # Build VTK unstructured Grid
#     if len(cell_types) == 1:
#         vtk_obj = _vtk.vtkUnstructuredGrid()
#         vtk_obj.SetPoints(vtkpoints)
#         vtk_obj.SetCells(cell_types[0], cell_array)
#     else:
#         cell_types_vtk = numpy_support.numpy_to_vtk(cell_types, deep=True, array_type=uchartype_vtk)
#         offsets_vtk = numpy_support.numpy_to_vtk(offsets, deep=True, array_type=idtype_vtk)
#         vtk_obj = _vtk.vtkUnstructuredGrid()
#         vtk_obj.SetPoints(vtkpoints)
#         vtk_obj.SetCells(cell_types_vtk, offsets_vtk, cell_array)

#     elemdata = {
#         "nodes": vtk_obj.GetPointData(),
#         "cells": vtk_obj.GetCellData(),
#     }
#     # And add the variables
#     for var in cellsvar:
#         np_array = cellsvar[var]
#         # cast without copy if possible
#         if np_array.flags.contiguous:
#             np_array = np_array.astype(np.float64, copy=False)
#         else:
#             np_array = np_array.astype(np.float64)
#         # zero-copy if contiguous array
#         vtkarray = numpy_support.numpy_to_vtk(np_array, deep=False)
#         vtkarray.charlesx_array = np_array  # avoid garbage collection
#         vtkarray.SetName(var)
#         #
#         elemdata["cells"].AddArray(vtkarray)
#     for var in nodesvar:
#         np_array = nodesvar[var]
#         # cast without copy if possible
#         if np_array.flags.contiguous:
#             np_array = np_array.astype(np.float64, copy=False)
#         else:
#             np_array = np_array.astype(np.float64)
#         # zero-copy if contiguous array
#         vtkarray = numpy_support.numpy_to_vtk(np_array, deep=False)
#         vtkarray.charlesx_array = np_array  # avoid garbage collection
#         vtkarray.SetName(var)
#         #
#         elemdata["nodes"].AddArray(vtkarray)
#     # Finally, write the file
#     writer = _vtk.vtkXMLUnstructuredGridWriter()
#     if _vtk.vtkVersion.GetVTKMajorVersion() >= 6:
#         writer.SetInputData(vtk_obj)
#     else:
#         writer.SetInput(vtk_obj)
#     filename = "restart.%d"%(simulation_state["step"])
#     writer.SetFileName(filename + '.vtu')
#     writer.SetDataModeToBinary()
#     writer.Write()
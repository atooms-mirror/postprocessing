import _neighbors_module
import f90wrap.runtime
import logging

class Neighbors_Module(f90wrap.runtime.FortranModule):
    """
    Module neighbors_module
    
    
    Defined at neighbors.f90 lines 1-43
    
    """
    @staticmethod
    def neighbors(box, pos, ids, rcut, nn, neigh):
        """
        neighbors(box, pos, ids, rcut, nn, neigh)
        
        
        Defined at neighbors.f90 lines 13-43
        
        Parameters
        ----------
        box : float array
        pos : float array
        ids : int array
        rcut : float array
        nn : int array
        neigh : int array
        
        """
        _neighbors_module.f90wrap_neighbors(box=box, pos=pos, ids=ids, rcut=rcut, nn=nn, \
            neigh=neigh)
    
    _dt_array_initialisers = []
    

neighbors_module = Neighbors_Module()


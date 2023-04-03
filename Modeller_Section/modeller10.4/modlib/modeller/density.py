"""Classes for handling EM (electron microscopy) density files."""

import _modeller
from modeller import modfile
from modeller.util.modobject import ModObject
from modeller.util import array
from modeller.util.deprecation import _deprecation_handler

__docformat__ = "epytext en"


class Density(ModObject):
    """Holds all information from an EM (electron microscopy) density file"""
    _modpt = None
    __read_vars = None
    __free_func = _modeller.mod_density_free
    env = None

    def __init__(self, env, **vars):
        self._modpt = _modeller.mod_density_new(self)
        self.env = env.copy()
        if len(vars) > 0:
            self.read(**vars)

    def __setstate__(self, d):
        self.__dict__.update(d)
        self._modpt = _modeller.mod_density_new(self)
        if self.__read_vars is not None:
            self.read(**self.__read_vars)

    def __del__(self):
        if self._modpt:
            self.__free_func(self._modpt)

    def read(self, file, resolution, em_map_size=0, voxel_size=0.,
             em_density_format='XPLOR', filter_type='NONE',
             filter_values=(0., 0.), density_type='SPHERE',
             px=None, py=None, pz=None, cc_func_type='CCF'):
        """Read in a density file"""
        self.__read_vars = {}
        for var in ('file', 'em_map_size', 'voxel_size', 'resolution',
                    'em_density_format', 'filter_type', 'filter_values',
                    'density_type', 'px', 'py', 'pz', 'cc_func_type'):
            self.__read_vars[var] = eval(var)
        read_origin = False
        if px is None or py is None or pz is None:
            px = py = pz = 0.
            read_origin = True
        fh = modfile._get_filehandle(file, 'rb')
        return _modeller.mod_density_read(
            self._modpt, fh.file_pointer, em_density_format, em_map_size,
            filter_type, voxel_size, resolution, filter_values, density_type,
            px, py, pz, read_origin, cc_func_type)

    def grid_search(self, em_pdb_name, chains_num, em_density_format='XPLOR',
                    num_structures=1, dock_order='INPUT', start_type='CENTER',
                    translate_type='NONE', number_of_steps=1,
                    angular_step_size=0, temperature=293.0,
                    best_docked_models=1, em_fit_output_file='em_fit.out'):
        """Dock a structure into the EM density map"""
        return _modeller.mod_density_grid_search(
            self._modpt, em_density_format, num_structures, dock_order,
            start_type, translate_type, number_of_steps, angular_step_size,
            temperature, best_docked_models, em_fit_output_file, em_pdb_name,
            chains_num)

    def __get_resolution(self):
        return _modeller.mod_density_resolution_get(self._modpt)

    def __set_resolution(self, val):
        _modeller.mod_density_resolution_set(self._modpt, val)

    def __get_sigma_factor(self):
        return _modeller.mod_density_sigma_factor_get(self._modpt)

    def __set_sigma_factor(self, val):
        _modeller.mod_density_sigma_factor_set(self._modpt, val)

    def __get_grid_sqr_sum(self):
        return _modeller.mod_density_grid_sqr_sum_get(self._modpt)

    def __get_grid(self):
        return array.Float3DArray(_modeller.mod_density_grid_get(self._modpt),
                                  (_modeller.mod_density_nz_get(self._modpt),
                                   _modeller.mod_density_ny_get(self._modpt),
                                   _modeller.mod_density_nx_get(self._modpt)))

    def __get_norm_factor(self):
        return _modeller.mod_density_norm_factor_get(self._modpt)

    def __get_voxel_size(self):
        return _modeller.mod_density_vox_size_get(self._modpt)

    def __get_px(self):
        return _modeller.mod_density_px_get(self._modpt)

    def __get_py(self):
        return _modeller.mod_density_py_get(self._modpt)

    def __get_pz(self):
        return _modeller.mod_density_pz_get(self._modpt)

    resolution = property(__get_resolution, __set_resolution, doc="Resolution")
    sigma_factor = property(__get_sigma_factor, __set_sigma_factor,
                            doc="Sigma factor")
    grid_sqr_sum = property(__get_grid_sqr_sum,
                            doc="Sum of squares over all grid points")
    grid = property(__get_grid, doc="3D grid of density values (z,y,x)")
    norm_factor = property(__get_norm_factor, doc="Normalization factor")
    voxel_size = property(__get_voxel_size, doc="Voxel size")
    px = property(__get_px, doc="x coordinate of origin")
    py = property(__get_py, doc="y coordinate of origin")
    pz = property(__get_pz, doc="z coordinate of origin")


# Modeller 9 compatibility
class density(Density):
    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(density)
        Density.__init__(self, *args, **keys)

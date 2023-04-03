"""Classes for handling arbitrary sets of atoms."""

import _modeller
import sys
from modeller import physical, GroupRestraints, EnergyData, modfile
from modeller.energy_profile import EnergyProfile
from modeller.util import modutil
from modeller.util.matrix import matrix_to_list
from modeller.util.deprecation import _deprecation_handler

__docformat__ = "epytext en"


class SelectionIterator(object):
    """Utility class to iterate over all atoms in a selection."""

    def __init__(self, mdl, seldict):
        self.mdl = mdl
        self.seliter = seldict.__iter__()

    def __iter__(self):
        return self

    if sys.version_info[0] >= 3:
        def __next__(self):
            obj = next(self.seliter)
            return self.mdl.atoms[obj - 1]
    else:
        def next(self):
            obj = self.seliter.next()
            return self.mdl.atoms[obj - 1]


class Selection(object):
    """An arbitrary set of atoms. 'atom' or 'residue' objects can be
       added to or removed from the selection, and selections can be
       manipulated in the same way as Python-standard sets -
       e.g. :meth:`union`, :meth:`intersection`, :meth:`difference`, or the
       equivalent ``|``, ``&`` and ``-`` operators.  Note that ``obj in sel``
       is only True if ALL atoms in ``obj`` are in the selection ``sel``.
       If you want to check for partial selection, use
       ``len(sel.intersection([obj])) > 0``."""

    # Global reference to the DOPE restraints files, so that they're loaded
    # only once
    dope_restraints = None
    dopehr_restraints = None

    def __init__(self, *atoms):
        self.__selection = {}
        self.__mdl = None
        self.union_update(atoms)

    def get_atom_indices(self):
        """Get the integer indices of all atoms in this selection

           :return: atom indices
           :rtype: list of ints
        """
        if sys.version_info[:2] == (2, 3):
            keys = self.__selection.keys()
            keys.sort()
        else:
            keys = sorted(self.__selection.keys())
        return (keys, self.__mdl)

    def get_model(self):
        """Get the model object which all selected atoms belong to.

           :return: the model object
           :rtype: :class:`Model`"""
        return self.__mdl

    def __repr__(self):
        suffix = "s"
        num = len(self.__selection)
        if num == 1:
            suffix = ""
        return "Selection of %d atom%s" % (num, suffix)

    def __str__(self):
        return "<%s>" % repr(self)

    def __len__(self):
        return len(self.__selection)

    def __check_object(self, obj):
        if not hasattr(obj, 'get_atom_indices'):
            raise TypeError(("Invalid object %s for selection - try atom, " +
                             "residue, model, Selection objects") % repr(obj))

    def __contains__(self, obj):
        self.__check_object(obj)
        (inds, mdl) = obj.get_atom_indices()
        if mdl is not self.__mdl:
            return False
        for ind in inds:
            if ind not in self.__selection:
                return False
        return True

    def __iter__(self):
        return SelectionIterator(self.__mdl, self.__selection)

    def __typecheck(self, obj):
        if not isinstance(obj, Selection):
            raise TypeError("Must use Selection objects here")

    def __typecoerce(self, obj):
        if isinstance(obj, Selection):
            return obj
        else:
            return Selection(*obj)

    def issubset(self, t):
        """Tests whether every element is also in t"""
        t = self.__typecoerce(t)
        for obj in self:
            if obj not in t:
                return False
        return True

    def __le__(self, t):
        self.__typecheck(t)
        return self.issubset(t)

    def issuperset(self, t):
        """Tests whether every element in t is also in this set"""
        t = self.__typecoerce(t)
        for obj in t:
            if obj not in self:
                return False
        return True

    def __ge__(self, t):
        self.__typecheck(t)
        return self.issuperset(t)

    def __eq__(self, t):
        self.__typecheck(t)
        return self.issubset(t) and t.issubset(self)

    def __lt__(self, t):
        self.__typecheck(t)
        return self.issubset(t) and not t.issubset(self)

    def __gt__(self, t):
        self.__typecheck(t)
        return t.issubset(self) and not self.issubset(t)

    def add(self, obj):
        """Adds a new object (e.g. atom, residue) to the selection"""
        if not hasattr(obj, "get_atom_indices") \
           and modutil.non_string_iterable(obj):
            for x in obj:
                self.add(x)
        else:
            self.__check_object(obj)
            (inds, mdl) = obj.get_atom_indices()
            if self.__mdl is None:
                self.__mdl = mdl
            elif self.__mdl is not mdl:
                raise ValueError("All atoms must be in the same model!")
            for ind in inds:
                self.__selection[ind] = None

    def remove(self, obj):
        """Removes an object (e.g. atom, residue) from the selection;
           raises KeyError if not present"""
        if not hasattr(obj, "get_atom_indices") \
           and modutil.non_string_iterable(obj):
            for x in obj:
                self.remove(x)
        else:
            if obj not in self:
                raise KeyError(obj)
            (inds, mdl) = obj.get_atom_indices()
            if self.__mdl is not mdl:
                raise ValueError("All atoms must be in the same model!")
            for ind in inds:
                del self.__selection[ind]

    def discard(self, obj):
        """Removes an object (e.g. atom, residue) from the selection,
           if present"""
        try:
            self.remove(obj)
        except KeyError:
            pass

    def union(self, t):
        """Get a new selection containing objects from both selections
           (the `|` operator is equivalent).

           :param t: the other selection to use for the union.
           :type t: :class:`Selection`
           :return: the new selection
           :rtype: :class:`Selection`
        """
        newobj = Selection()
        # Handle the case where both selections are empty
        newobj.__mdl = self.__mdl
        for obj in self:
            newobj.add(obj)
        for obj in t:
            newobj.add(obj)
        return newobj

    def __or__(self, t):
        self.__typecheck(t)
        return self.union(t)

    def intersection(self, t):
        """Create and return a new selection containing objects common to both
           selections (the `&` operator is equivalent).

           :param t: the other selection to use for the intersection
           :type t: :class:`Selection`
           :return: the new selection
           :rtype: :class:`Selection`"""
        t = self.__typecoerce(t)
        newobj = Selection()
        # Handle the case where both selections are empty
        newobj.__mdl = self.__mdl
        for obj in self:
            if obj in t:
                newobj.add(obj)
        return newobj

    def __and__(self, t):
        self.__typecheck(t)
        return self.intersection(t)

    def difference(self, t):
        """Create and return a new selection containing objects in this
           selection but not in ``t``  (the `-` operator is equivalent).

           :param t: the other selection to use for the difference
           :type t: :class:`Selection`
           :return: the new selection
           :rtype: :class:`Selection`"""
        t = self.__typecoerce(t)
        newobj = Selection()
        # Handle the case where both selections are empty
        newobj.__mdl = self.__mdl
        for obj in self:
            if obj not in t:
                newobj.add(obj)
        return newobj

    def __sub__(self, t):
        self.__typecheck(t)
        return self.difference(t)

    def symmetric_difference(self, t):
        """Return a new selection containing objects in either this selection
           or ``t``, but not both (the `^` operator is equivalent).

           :param t: the other selection to use for the symmetric difference
           :type t: :class:`Selection`
           :return: the new selection
           :rtype: :class:`Selection`"""
        t = self.__typecoerce(t)
        newobj = Selection()
        # Handle the case where both selections are empty
        newobj.__mdl = self.__mdl
        for obj in self:
            if obj not in t:
                newobj.add(obj)
        for obj in t:
            if obj not in self:
                newobj.add(obj)
        return newobj

    def __xor__(self, t):
        self.__typecheck(t)
        return self.symmetric_difference(t)

    def copy(self):
        """Make and return a copy of this selection.

           :return: a copy of this selection.
           :rtype: :class:`Selection`"""
        newobj = Selection()
        newobj.__selection = self.__selection.copy()
        newobj.__mdl = self.__mdl
        return newobj

    def union_update(self, t):
        """Like :meth:`union`, but this selection is updated in place with
           the result (the `|=` operator is equivalent)."""
        for obj in t:
            self.add(obj)
        return self

    def __ior__(self, t):
        self.__typecheck(t)
        return self.union_update(t)

    def intersection_update(self, t):
        """Like :meth:`intersection`, but this selection is updated in
           place with the result (the `&=` operator is equivalent)."""
        t = self.__typecoerce(t)
        c = self.copy()
        for obj in c:
            if obj not in t:
                self.remove(obj)
        return self

    def __iand__(self, t):
        self.__typecheck(t)
        return self.intersection_update(t)

    def difference_update(self, t):
        """Like :meth:`difference`, but this selection is updated in
           place with the result (the `-=` operator is equivalent)."""
        t = self.__typecoerce(t)
        c = self.copy()
        for obj in c:
            if obj in t:
                self.remove(obj)
        return self

    def __isub__(self, t):
        self.__typecheck(t)
        return self.difference_update(t)

    def symmetric_difference_update(self, t):
        """Like :meth:`symmetric_difference`, but this selection is updated in
           place with the result (the `^=` operator is equivalent)."""
        t = self.__typecoerce(t)
        c = self.copy()
        for obj in t:
            if obj in c:
                self.remove(obj)
            else:
                self.add(obj)
        return self

    def __ixor__(self, t):
        self.__typecheck(t)
        return self.symmetric_difference_update(t)

    def clear(self):
        """Remove all atoms from the selection"""
        self.__mdl = None
        self.__selection = {}

    def __require_indices(self):
        """Get atom indices, and fail if there are none"""
        (inds, mdl) = self.get_atom_indices()
        if mdl is None or len(inds) == 0:
            raise ValueError("Selection contains no atoms")
        return (inds, mdl)

    def translate(self, vector):
        """Translate the selection by the given vector"""
        (inds, mdl) = self.__require_indices()
        _modeller.mod_selection_translate(mdl.modpt, inds, vector)

    def rotate_origin(self, axis, angle):
        """Rotate the selection about the given axis through the origin,
           by the given angle (in degrees)"""
        (inds, mdl) = self.__require_indices()
        if sum([x * x for x in axis]) < 1e-8:
            raise ValueError("axis must be a vector of non-zero length")
        _modeller.mod_selection_rotate_axis(mdl.modpt, inds, axis, angle)

    def rotate_mass_center(self, axis, angle):
        """Rotate the selection about the given axis through the mass center,
           by the given angle (in degrees)"""
        com = self.mass_center
        self.translate([-a for a in com])
        try:
            self.rotate_origin(axis, angle)
        finally:
            # If rotate_origin fails, at least restore the
            # original conformation
            self.translate(com)

    def transform(self, matrix):
        """Transform the selection coordinates with the given matrix"""
        (inds, mdl) = self.__require_indices()
        lst = matrix_to_list(matrix)
        _modeller.mod_selection_transform(mdl.modpt, inds, lst)

    def find_atoms(self, restyp, atom_names, min_selected):
        if min_selected <= 0:
            raise ValueError("min_selected must be at least 1")
        mdl = self.get_model()
        if mdl is None or len(self) == 0:
            raise ValueError("Selection contains no atoms")
        from modeller import model_topology
        return model_topology.FindAtoms(mdl, restyp, atom_names,
                                        self.__check_selected, min_selected)

    def find_chi1_dihedrals(self, restyp, min_selected):
        return self.__find_dihedrals(restyp, 5, min_selected)

    def find_chi2_dihedrals(self, restyp, min_selected):
        return self.__find_dihedrals(restyp, 6, min_selected)

    def find_chi3_dihedrals(self, restyp, min_selected):
        return self.__find_dihedrals(restyp, 7, min_selected)

    def find_chi4_dihedrals(self, restyp, min_selected):
        return self.__find_dihedrals(restyp, 8, min_selected)

    def __find_dihedrals(self, restyp, dihedral_type, min_selected):
        if min_selected <= 0:
            raise ValueError("min_selected must be at least 1")
        mdl = self.get_model()
        if mdl is None or len(self) == 0:
            raise ValueError("Selection contains no atoms")
        from modeller import model_topology
        return model_topology.FindDihedrals(mdl, restyp, dihedral_type,
                                            self.__check_selected,
                                            min_selected)

    def __check_selected(self, atoms, min_selected):
        num_selected = 0
        min_selected = min(min_selected, len(atoms))
        for a in atoms:
            if a in self:
                num_selected += 1
                if num_selected >= min_selected:
                    return True
        return False

    def write(self, file, model_format='PDB', no_ter=False, extra_data=""):
        """Write selection coordinates to a file"""
        (inds, mdl) = self.__require_indices()
        fh = modfile._get_filehandle(file, 'w')
        return _modeller.mod_model_write(mdl.modpt, mdl.env.libs.modpt, inds,
                                         fh.file_pointer, model_format,
                                         no_ter, False, extra_data)

    def by_residue(self):
        """Create and return a new selection, in which any residues in the
           existing selection that have at least one selected atom are now
           entirely selected.

           :return: the new selection
           :rtype: :class:`selection`"""
        return self.extend_by_residue(0)

    def extend_by_residue(self, extension):
        """Create and return a new selection, in which any residues with
           at least one selected atom in the existing selection are now
           entirely selected. Additionally, ``extension`` residues around
           each selected residue are selected.

           :param int extension: the number of residues to extend the
                  selection by
           :return: the new selection
           :rtype: :class:`Selection`"""
        newobj = self.copy()
        (inds, mdl) = self.get_atom_indices()
        if mdl is not None:
            for (n, res) in enumerate(mdl.residues):
                if len(self.intersection([res])) > 0:
                    for i in range(max(0, n-extension),
                                   min(len(mdl.residues), n+extension+1)):
                        newobj.add(mdl.residues[i])
        return newobj

    def only_sidechain(self):
        """Create and return a new selection, containing only sidechain
           atoms from the current selection.

           :return: the new selection containing only sidechain atoms
           :rtype: :class:`Selection`"""
        return self - self.only_mainchain()

    def only_mainchain(self):
        """Create and return a new selection, containing only mainchain
           atoms from the current selection.

           :return: the new selection containing only mainchain atoms
           :rtype: :class:`Selection`"""
        return self.only_atom_types('O OT1 OT2 OXT C CA N')

    def only_atom_types(self, atom_types):
        """Create and return a new selection, containing only atoms from
           the current selection of the given type(s).

           :param str atom_types: A space-separated list of type(s)
                  (e.g. 'CA CB').
           :return: A new selection containing only atoms of the given type(s).
           :rtype: :class:`Selection`"""
        if modutil.non_string_iterable(atom_types):
            atom_types = " ".join(atom_types)
        return self.__filter(_modeller.mod_selection_atom_types, atom_types)

    def only_residue_types(self, residue_types):
        """Create and return a new selection, containing only atoms from
           the current selection in residues of the given type(s).

           :param str residue_types: A space-separated list of type(s)
                  (e.g. 'ALA ASP').
           :return: A new selection containing only atoms in residues
                  of the given type(s).
           :rtype: :class:`Selection`"""
        if modutil.non_string_iterable(residue_types):
            residue_types = " ".join(residue_types)
        libs = self.__get_libs()
        return self.__filter(_modeller.mod_selection_residue_types,
                             residue_types, libs)

    def only_std_residues(self):
        """Create and return a new selection, containing only atoms in
           standard residue types (i.e. everything but HETATM).

           :return: The new selection
           :rtype: :class:`Selection`"""
        libs = self.__get_libs()
        return self.__filter(_modeller.mod_selection_std_residues, libs)

    def only_no_topology(self):
        """Create and return a new selection, containing only atoms in
           residues that have no defined topology.

           :return: The new selection
           :rtype: :class:`Selection`"""
        libs = self.__get_libs()
        return self.__filter(_modeller.mod_selection_no_topology, libs)

    def only_het_residues(self):
        """Create and return a new selection, containing only atoms in
           HETATM residues.

           :return: The new selection
           :rtype: :class:`Selection`"""
        libs = self.__get_libs()
        return self.__filter(_modeller.mod_selection_het_residues, libs)

    def only_water_residues(self):
        """Create and return a new selection, containing only atoms in
           water residues.

           :return: The new selection
           :rtype: :class:`Selection`"""
        libs = self.__get_libs()
        return self.__filter(_modeller.mod_selection_water_residues, libs)

    def only_defined(self):
        """Create and return a new selection, containing only atoms with defined
           coordinates.

           :return: The new selection
           :rtype: :class:`Selection`"""
        return self.__filter(_modeller.mod_selection_defined)

    def select_sphere(self, radius):
        """Create and return a new selection, containing all atoms within
           the given distance from any atom in the current selection.

           :param float radius: Distance in angstroms
           :return: The new selection
           :rtype: :class:`Selection`"""
        s = Selection()
        for x in self:
            s.add(x.select_sphere(radius))
        return s

    def mutate(self, residue_type):
        """Mutate selected residues"""
        (inds, mdl) = self.__require_indices()
        return _modeller.mod_selection_mutate(mdl.modpt, mdl.env.libs.modpt,
                                              inds, residue_type)

    def randomize_xyz(self, deviation):
        """Randomize coordinates"""
        (inds, mdl) = self.__require_indices()
        return _modeller.mod_randomize_xyz(mdl.modpt, mdl.env.libs.modpt, inds,
                                           deviation)

    def superpose(self, mdl2, aln, fit=True, superpose_refine=False,
                  rms_cutoff=3.5, reference_atom='', reference_distance=3.5,
                  refine_local=True, swap_atoms_in_res=''):
        """Superpose the input model on this selection, given an alignment of
           the models"""
        from modeller import superpose
        (inds, mdl) = self.__require_indices()
        retval = _modeller.mod_superpose(mdl.modpt, mdl2.modpt, aln.modpt,
                                         mdl.env.libs.modpt, inds,
                                         swap_atoms_in_res, reference_atom,
                                         reference_distance, superpose_refine,
                                         fit, refine_local, rms_cutoff)
        return superpose.SuperposeData(*retval)

    def rotate_dihedrals(self, deviation, change,
                         dihedrals=('PHI', 'PSI', 'CHI1', 'CHI2',
                                    'CHI3', 'CHI4')):
        """Optimize or randomize dihedral angles"""
        (inds, mdl) = self.__require_indices()
        return _modeller.mod_rotate_dihedrals(mdl.modpt, mdl.env.libs.modpt,
                                              inds, deviation, change,
                                              dihedrals)

    def hot_atoms(self, pick_hot_cutoff, residue_span_range=(0, 99999),
                  viol_report_cut=physical.Values(
                      default=4.5, chi1_dihedral=999, chi2_dihedral=999,
                      chi3_dihedral=999, chi4_dihedral=999,
                      chi5_dihedral=999, phi_psi_dihedral=6.5,
                      nonbond_spline=999, accessibility=999,
                      density=999, gbsa=999, em_density=999),
                  schedule_scale=None, edat=None):
        """Return a new selection containing all atoms violating restraints.

           :return: The new selection
           :rtype: :class:`Selection`"""
        (inds, mdl) = self.__require_indices()

        if edat is None:
            edat = mdl.env.edat
        if schedule_scale is None:
            schedule_scale = mdl.env.schedule_scale
        func = _modeller.mod_selection_hot_atoms
        newinds = func(mdl.modpt, edat.modpt, mdl.env.libs.modpt, inds,
                       residue_span_range, pick_hot_cutoff, viol_report_cut,
                       schedule_scale)
        newobj = Selection()
        newobj.__mdl = mdl
        newobj.__selection = dict.fromkeys(newinds)
        return newobj

    def objfunc(self, edat=None, residue_span_range=(0, 99999),
                schedule_scale=physical.Values(default=1.0)):
        """Get just the objective function value, without derivatives"""
        (inds, mdl) = self.__require_indices()
        if edat is None:
            edat = mdl.env.edat
        return _modeller.mod_selection_objfunc(mdl.modpt, edat.modpt,
                                               mdl.env.libs.modpt, inds,
                                               residue_span_range,
                                               schedule_scale)

    def energy(self, asgl_output=False, normalize_profile=False,
               residue_span_range=(0, 99999), output='LONG', file='default',
               viol_report_cut=physical.Values(
                   default=4.5, chi1_dihedral=999, chi2_dihedral=999,
                   chi3_dihedral=999, chi4_dihedral=999, chi5_dihedral=999,
                   phi_psi_dihedral=6.5, nonbond_spline=999, accessibility=999,
                   density=999, gbsa=999, em_density=999),
               viol_report_cut2=physical.Values(default=2.0),
               smoothing_window=3, schedule_scale=None, edat=None):
        """Evaluate the objective function given restraints"""
        (inds, mdl) = self.__require_indices()

        if edat is None:
            edat = mdl.env.edat
        if schedule_scale is None:
            schedule_scale = mdl.env.schedule_scale
        func = _modeller.mod_energy
        (molpdf, terms) = func(mdl.modpt, edat.modpt, mdl.env.libs.modpt, inds,
                               asgl_output, normalize_profile,
                               residue_span_range, output, file,
                               smoothing_window, viol_report_cut,
                               viol_report_cut2, schedule_scale)
        terms = physical.from_list(terms)
        return (molpdf, terms)

    def get_dope_potential(self):
        """Get the spline data for the DOPE statistical potential"""
        (inds, mdl) = self.__require_indices()
        if not Selection.dope_restraints:
            Selection.dope_restraints = \
                GroupRestraints(mdl.env, classes='${LIB}/atmcls-mf.lib',
                                parameters='${LIB}/dist-mf.lib')
        return Selection.dope_restraints

    def get_dopehr_potential(self):
        """Get the spline data for the DOPE-HR statistical potential"""
        (inds, mdl) = self.__require_indices()
        if not Selection.dopehr_restraints:
            Selection.dopehr_restraints = \
                GroupRestraints(mdl.env, classes='${LIB}/atmcls-mf.lib',
                                parameters='${LIB}/dist-mfhr.lib')
        return Selection.dopehr_restraints

    def get_dope_energy_data(self):
        """Get ideal energy_data terms for DOPE potential evaluations"""
        return EnergyData(contact_shell=15.0, dynamic_modeller=True,
                          dynamic_lennard=False, dynamic_sphere=False,
                          excl_local=(False, False, False, False))

    def assess_dope(self, **vars):
        """Assess the selection with the DOPE potential"""
        return self._dope_energy(self.get_dope_potential(), "DOPE", **vars)

    def assess_dopehr(self, **vars):
        """Assess the selection with the DOPE-HR potential"""
        return self._dope_energy(
            self.get_dopehr_potential(), "DOPE-HR", **vars)

    def get_energy_profile(self, physical_type, edat=None):
        """Get a per-residue energy profile, plus the number of restraints on
           each residue, and the RMS minimum and heavy violations"""
        (inds, mdl) = self.__require_indices()
        if edat is None:
            edat = mdl.env.edat
        scaln = physical.Values(default=0.)
        scaln[physical_type] = 1.
        prof = _modeller.mod_rms_profile(mdl.modpt, edat.modpt,
                                         mdl.env.libs.modpt, inds, (1, 9999),
                                         True, False, physical_type.get_type(),
                                         scaln)
        return EnergyProfile(*prof)

    def get_dope_profile(self):
        """Get a per-residue DOPE energy profile"""
        return self._dope_profile(self.get_dope_potential())

    def get_dopehr_profile(self):
        """Get a per-residue DOPE-HR energy profile"""
        return self._dope_profile(self.get_dopehr_potential())

    def _dope_profile(self, gprsr):
        """Internal function to get a DOPE or DOPE-HR energy profile"""
        mdl = self.__mdl
        edat = self.get_dope_energy_data()
        old_gprsr = mdl.group_restraints
        mdl.group_restraints = gprsr
        try:
            prof = self.get_energy_profile(physical.nonbond_spline, edat)
        finally:
            mdl.group_restraints = old_gprsr
        return prof

    def assess(self, assessor, output='SHORT NO_REPORT', **vars):
        """Assess with the given assessor object
           (e.g. :class:`soap_loop.Scorer`)."""
        print(">> Model assessment by %s" % assessor.name)
        molpdf, terms = assessor._assess(self, output=output, **vars)
        print("%s                     : %12.6f" % (assessor.name, molpdf))
        return molpdf

    def _dope_energy(self, gprsr, name, output='SHORT NO_REPORT',
                     residue_span_range=(1, 9999),
                     schedule_scale=physical.Values(
                         default=0., nonbond_spline=1.), **vars):
        """Internal function to do DOPE or DOPE-HR assessment"""
        mdl = self.__mdl
        print(">> Model assessment by %s potential" % name)
        edat = self.get_dope_energy_data()
        old_gprsr = mdl.group_restraints
        mdl.group_restraints = gprsr
        try:
            (molpdf, terms) = \
                self.energy(edat=edat, residue_span_range=residue_span_range,
                            output=output, schedule_scale=schedule_scale,
                            **vars)
        finally:
            mdl.group_restraints = old_gprsr
        print("%s score               : %12.6f" % (name, molpdf))
        return molpdf

    def debug_function(self, residue_span_range=(0, 99999),
                       debug_function_cutoff=(0.01, 0.001, 0.1),
                       detailed_debugging=False, schedule_scale=None,
                       edat=None):
        """Test code self-consistency"""
        (inds, mdl) = self.__require_indices()

        if edat is None:
            edat = mdl.env.edat
        if schedule_scale is None:
            schedule_scale = mdl.env.schedule_scale
        func = _modeller.mod_debug_function
        return func(mdl.modpt, edat.modpt, mdl.env.libs.modpt, inds,
                    residue_span_range, debug_function_cutoff,
                    detailed_debugging, schedule_scale)

    def unbuild(self):
        """Undefine all coordinates"""
        (inds, mdl) = self.__require_indices()
        return _modeller.mod_selection_unbuild(mdl=mdl.modpt, sel1=inds)

    def __get_libs(self):
        if self.__mdl:
            return self.__mdl.env.libs.modpt
        else:
            return None

    def __filter(self, func, *args):
        newobj = Selection()
        (inds, mdl) = self.get_atom_indices()
        if mdl is not None:
            newinds = func(mdl.modpt, inds, *args)
            newobj.__mdl = mdl
            newobj.__selection = dict.fromkeys(newinds)
        return newobj

    def __get_com(self):
        (inds, mdl) = self.__require_indices()
        return _modeller.mod_selection_com_get(mdl.modpt, inds)

    def __set_com(self, val):
        com = self.mass_center
        self.translate([val[0] - com[0], val[1] - com[1], val[2] - com[2]])

    def __get_x(self):
        return self.__get_com()[0]

    def __get_y(self):
        return self.__get_com()[1]

    def __get_z(self):
        return self.__get_com()[2]

    def __set_x(self, val):
        self.translate([val - self.x, 0, 0])

    def __set_y(self, val):
        self.translate([0, val - self.y, 0])

    def __set_z(self, val):
        self.translate([0, 0, val - self.z])

    mass_center = property(__get_com, __set_com,
                           doc="Coordinates of mass center")
    x = property(__get_x, __set_x, doc="x coordinate of mass center")
    y = property(__get_y, __set_y, doc="y coordinate of mass center")
    z = property(__get_z, __set_z, doc="z coordinate of mass center")


# Modeller 9 compatibility
class selection(Selection):
    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(selection)
        Selection.__init__(self, *args, **keys)

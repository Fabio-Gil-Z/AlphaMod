"""Classes for handling user-defined energy terms"""

from modeller.util.modlist import LinkList
from modeller.util.deprecation import _deprecation_handler
import _modeller


class EnergyTerm(object):
    """Base class for user-defined energy terms. Override the :meth:`eval`
       method and the :attr:`_physical_type` member in your own subclasses."""

    _physical_type = None
    _cutoff = 0.0

    def _add_term(self, edat, indx):
        """Register this energy term with Modeller"""
        _modeller.mod_energy_term_new(edat, indx, self.eval, self._cutoff,
                                      self._physical_type.get_type())

    def eval(self, mdl, deriv, indats):
        raise NotImplementedError("No Python function for energy term")

    def indices_to_atoms(self, mdl, atom_indices):
        """Converts Modeller-style atom indices into atom objects"""
        return [mdl.atoms[x-1] for x in atom_indices]


class TermList(LinkList):
    """A list of :class:`energy_term` objects"""
    def __init__(self, edat):
        self.__edat = edat
        self.__list = []
        LinkList.__init__(self)

    def __setstate__(self, d):
        self.__dict__.update(d)
        for (indx, obj) in enumerate(self.__list):
            obj._add_term(self.__edat(), indx)

    def _insfunc(self, indx, obj):
        obj._add_term(self.__edat(), indx)
        self.__list.insert(indx, obj)

    def __len__(self):
        return len(self.__list)

    def _getfunc(self, indx):
        return self.__list[indx]

    def _delfunc(self, indx):
        del self.__list[indx]
        _modeller.mod_energy_term_del(self.__edat(), indx)


class AssessEnergyTerm(EnergyTerm):
    """Base class for an energy term that can also be used for assessment."""

    _edat = None

    def _assess(self, atmsel, schedule_scale=None, **vars):
        if schedule_scale is None:
            schedule_scale = self._get_schedule_scale()
        return atmsel.energy(edat=self._get_energy_data_cached(),
                             schedule_scale=schedule_scale, **vars)

    def __call__(self, atmsel):
        return (self.name, atmsel.assess(self))

    def _get_energy_data(self):
        from modeller import EnergyData
        return EnergyData(contact_shell=-999, dynamic_sphere=False)

    def _get_energy_data_cached(self):
        if self._edat is None:
            self._edat = self._get_energy_data_all()
        return self._edat

    def _get_energy_data_all(self):
        edat = self._get_energy_data()
        edat.energy_terms.append(self)
        return edat

    def _get_schedule_scale(self):
        from modeller import physical
        schedule_scale = physical.Values(default=0.)
        schedule_scale[self._group] = 1.
        return schedule_scale


# Modeller 9 compatibility
class energy_term(EnergyTerm):
    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(energy_term)
        EnergyTerm.__init__(self, *args, **keys)

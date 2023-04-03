"""Classes for scoring models with GB/SA implicit solvation"""

from modeller import terms, physical
import _modeller

__docformat__ = "epytext en"


class Scorer(terms.EnergyTerm):
    """Score the model using GB/SA implicit solvation."""

    _physical_type = physical.gbsa

    def __init__(self, library='${LIB}/solv.lib', solvation_model=1,
                 cutoff=8.0):
        terms.EnergyTerm.__init__(self)
        self.__library = library
        self.__solvation_model = solvation_model
        self.__cutoff = cutoff

    def _add_term(self, edat, indx):
        _modeller.mod_gbsa_create(edat, indx, self._physical_type.get_type(),
                                  self.__library, self.__solvation_model,
                                  self.__cutoff)

"""Classes for scoring models with the SOAP-Protein-OD potential"""

from modeller import terms, physical
import _modeller

__docformat__ = "epytext en"


class Scorer(terms.AssessEnergyTerm):
    """Score the model using the SOAP-Protein-OD potential."""

    name = "SOAP-Protein-OD score"

    def __init__(self, library='${LIB}/soap_protein_od.hdf5',
                 group=physical.xy_distance):
        terms.EnergyTerm.__init__(self)
        self._group = group
        self.__library = library

    def _add_term(self, edat, indx):
        _modeller.mod_soap_od_create(edat, indx, self._group.get_type(),
                                     self.__library)

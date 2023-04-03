"""Classes for scoring models with the SOAP-Loop potential"""

from modeller import terms, physical
import _modeller

__docformat__ = "epytext en"


class Scorer(terms.AssessEnergyTerm):
    """Score the model using the SOAP-Loop potential."""

    name = "SOAP-Loop score"

    def __init__(self, library='${LIB}/soap_loop.hdf5',
                 group=physical.xy_distance):
        terms.EnergyTerm.__init__(self)
        self._group = group
        self.__library = library

    def _add_term(self, edat, indx):
        _modeller.mod_soap_od_create(edat, indx, self._group.get_type(),
                                     self.__library)

    def _get_energy_data(self):
        from modeller import EnergyData
        return EnergyData(contact_shell=-999, dynamic_sphere=False,
                          excl_local=[True, False, False, True])

"""Classes for scoring models with the SOAP-PP potential"""

from modeller import terms, physical
import _modeller

__docformat__ = "epytext en"


class PairScorer(terms.AssessEnergyTerm):
    """Score the model using the SOAP-PP pairwise term."""
    name = 'SOAP-PP pair score'

    def __init__(self, library='${LIB}/soap_pp_pair.hdf5',
                 group=physical.xy_distance):
        terms.EnergyTerm.__init__(self)
        self.__library = library
        self._group = group

    def _add_term(self, edat, indx):
        _modeller.mod_soap_pair_create(edat, indx, self._group.get_type(),
                                       self.__library)


class AtomScorer(terms.AssessEnergyTerm):
    """Score the model using the SOAP-PP atomic (accessibility) term."""
    name = 'SOAP-PP access score'

    def __init__(self, library='${LIB}/soap_pp_atom.hdf5',
                 group=physical.accessibility):
        terms.EnergyTerm.__init__(self)
        self.__library = library
        self._group = group

    def _add_term(self, edat, indx):
        _modeller.mod_soap_access_create(edat, indx, self._group.get_type(),
                                         self.__library)


class Assessor(terms.AssessEnergyTerm):
    name = 'SOAP-PP score'

    def __init__(self, pair_library='${LIB}/soap_pp_pair.hdf5',
                 pair_group=physical.xy_distance,
                 atom_library='${LIB}/soap_pp_atom.hdf5',
                 atom_group=physical.accessibility):
        self.pair_scorer = PairScorer(pair_library, pair_group)
        self.atom_scorer = AtomScorer(atom_library, atom_group)

    def _get_energy_data_all(self):
        edat = self._get_energy_data()
        edat.energy_terms.append(self.pair_scorer)
        edat.energy_terms.append(self.atom_scorer)
        return edat

    def _get_schedule_scale(self):
        from modeller import physical
        schedule_scale = physical.Values(default=0.)
        schedule_scale[self.pair_scorer._group] = 1.
        schedule_scale[self.atom_scorer._group] = 1.
        return schedule_scale

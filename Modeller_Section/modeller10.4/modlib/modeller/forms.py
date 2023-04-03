import _modeller
from modeller import physical
from modeller.util import modutil
from modeller.util.deprecation import _deprecation_handler


class RestraintForm(object):
    """Base class for all restraint forms"""

    _builtin_index = None
    _user_index = None
    _use_array = False
    ndim = None

    def __init__(self, group, features, modality, parameters):
        if '_builtin_index' not in self.__class__.__dict__ \
           and self._user_index is None:
            self.__class__._user_index = \
                _modeller.mod_user_form_new2(self.eval, self.vmin, self.vheavy,
                                             self.rvmin, self.rvheavy,
                                             self.min_mean, self.heavy_mean,
                                             self.get_range)
        if not isinstance(group, physical.PhysicalType):
            raise TypeError("group should be a PhysicalType object")
        self._group = group
        self._features = features
        self._modality = modality
        self._parameters = parameters

    def __get_restraint(self, mdl):
        feats = self._features
        if not isinstance(feats, (tuple, list)):
            feats = (feats,)
        modal = self._modality
        if not isinstance(modal, (tuple, list)):
            modal = (modal,)
        if self.ndim is not None and \
           (self.ndim != len(modal) or self.ndim != len(feats)):
            raise ValueError("Incorrect number of dimensions - did you "
                             "forget to call add_dimension() ?")
        iftyp = []
        atoms = []
        natoms = []
        for ft in feats:
            iftyp.append(ft.get_type())
            (at, ftmdl) = ft.get_atom_indices()
            if ftmdl is None:
                raise ValueError("Feature contains no atoms!")
            elif ftmdl != mdl:
                raise ValueError("Feature is defined on a different model!")
            natoms.append(len(at))
            atoms.extend(at)
        if isinstance(self._use_array, bool):
            prm = self._parameters
            use = self._use_array
        else:
            prm = (self._use_array,)
            use = False
        return (self.get_type(), self._group._num, iftyp, modal, natoms, atoms,
                prm, use)

    def _add_restraint(self, rsr, mdl):
        rsrdata = self.__get_restraint(mdl)
        ret = _modeller.mod_restraints_add(rsr.modpt, *rsrdata)
        if ret > 0:
            return ret

    def get_type(cls):
        """Get the numeric type of this form"""
        if cls._user_index is not None:
            return cls._user_index
        elif cls._builtin_index is not None:
            return cls._builtin_index
        else:
            raise ValueError("Cannot get type - no objects created yet")
    get_type = classmethod(get_type)

    def deltaf(self, feat, mean, iftyp):
        return _modeller.mod_feature_delta(feat, mean, iftyp)

    def eval(self, feats, iftyp, modal, deriv, param):
        raise NotImplementedError("Value not defined for form " + str(self))

    def vmin(self, feats, iftyp, modal, param):
        raise NotImplementedError("Min. violation not defined for form "
                                  + str(self))

    def vheavy(self, feats, iftyp, modal, param):
        raise NotImplementedError("Heavy violation not defined for form "
                                  + str(self))

    def rvmin(self, feats, iftyp, modal, param):
        raise NotImplementedError("Relative min. violation not defined " +
                                  "for form " + str(self))

    def rvheavy(self, feats, iftyp, modal, param):
        raise NotImplementedError("Relative heavy violation not defined " +
                                  "for form " + str(self))

    def min_mean(self, feats, iftyp, modal, param):
        raise NotImplementedError(
            "Min. mean not defined for form " + str(self))

    def heavy_mean(self, feats, iftyp, modal, param):
        raise NotImplementedError("Heavy mean not defined for form "
                                  + str(self))

    def get_range(self, iftyp, modal, param, spline_range):
        raise NotImplementedError("Feature range not defined for form "
                                  + str(self))


class LowerBound(RestraintForm):
    """Left Gaussian (harmonic lower bound)"""
    _builtin_index = 1

    def __init__(self, group, feature, mean, stdev):
        RestraintForm.__init__(self, group, feature, 0, (mean, stdev))


class UpperBound(RestraintForm):
    """Right Gaussian (harmonic upper bound)"""
    _builtin_index = 2

    def __init__(self, group, feature, mean, stdev):
        RestraintForm.__init__(self, group, feature, 0, (mean, stdev))


class Gaussian(RestraintForm):
    """Single Gaussian (harmonic potential)"""
    _builtin_index = 3

    def __init__(self, group, feature, mean, stdev):
        RestraintForm.__init__(self, group, feature, 0, (mean, stdev))


class MultiGaussian(RestraintForm):
    """Multiple Gaussian"""
    _builtin_index = 4

    def __init__(self, group, feature, weights, means, stdevs):
        lv = -1
        for var in (weights, means, stdevs):
            if (lv >= 0 and lv != len(var)) \
               or not isinstance(var, (tuple, list)):
                raise TypeError("weights, means and stdevs should all be "
                                + "sequences of the same length")
            lv = len(var)
        RestraintForm.__init__(self, group, feature, len(weights),
                               tuple(weights) + tuple(means) + tuple(stdevs))


class LennardJones(RestraintForm):
    """Lennard-Jones potential"""
    _builtin_index = 5

    def __init__(self, group, feature, A, B):
        RestraintForm.__init__(self, group, feature, 0, (A, B))


class Coulomb(RestraintForm):
    """Coulomb potential"""
    _builtin_index = 6

    def __init__(self, group, feature, q1, q2):
        RestraintForm.__init__(self, group, feature, 0, (q1, q2))


class Cosine(RestraintForm):
    """Cosine potential"""
    _builtin_index = 7

    def __init__(self, group, feature, phase, force, period):
        RestraintForm.__init__(self, group, feature, period, (phase, force))


class Factor(RestraintForm):
    """Simple scaling of feature value"""
    _builtin_index = 8

    def __init__(self, group, feature, factor):
        RestraintForm.__init__(self, group, feature, 0, (factor,))


class MultiBinormal(RestraintForm):
    """Multiple binormal"""
    _builtin_index = 9

    def __init__(self, group, features, weights, means, stdevs, correls):
        lv = -1
        for var in (weights, means, stdevs, correls):
            if (lv >= 0 and lv != len(var)) \
               or not isinstance(var, (tuple, list)):
                raise TypeError("weights, means, stdevs and correls should "
                                + "all be sequences of the same length")
            lv = len(var)
        for var in (means, stdevs):
            for elem in var:
                if not isinstance(elem, (tuple, list)) or len(elem) != 2:
                    raise TypeError("Each element of means and stdevs should "
                                    + "be a list of two elements")
        if not isinstance(features, (tuple, list)) or len(features) != 2:
            raise TypeError("features should be a list of two features")
        params = []
        for m in means:
            params += list(m)
        for m in stdevs:
            params += list(m)
        RestraintForm.__init__(self, group, features, [len(weights)]*2,
                               tuple(weights) + tuple(params) + tuple(correls))


class Spline(RestraintForm):
    """Cubic spline"""
    _builtin_index = 10
    def __init__(self, group, feature, open, low, high, delta, lowderiv,
                 highderiv, values, use_array=False):
        self._use_array = use_array
        if not isinstance(values, (tuple, list)):
            raise TypeError("values must be a sequence of spline points")
        spltyp = {False: -1.0, True: 0.0}
        RestraintForm.__init__(self, group, feature, len(values),
                               (spltyp[open], low, high, delta, lowderiv,
                                highderiv) + tuple(values))


class NDSpline(RestraintForm):
    """Multi-dimensional cubic spline"""

    _builtin_index = 10

    def __init__(self, group, values, dimensions=None, use_array=False):
        self._use_array = use_array
        modal = []
        param = []
        if dimensions is not None:
            modal = dimensions[:]
            nparam = 1
            for dim in modal:
                nparam *= dim
            if nparam != len(values):
                raise ValueError("Number of spline parameters is not equal" +
                                 " to the product of supplied dimensions")
            param = list(values)
        else:
            self._parse_values(values, modal, param, 0)
        RestraintForm.__init__(self, group, [], modal, param)
        self.ndim = 0

    def add_dimension(self, feature, open, low, high, delta, lowderiv,
                      highderiv):
        """Set up the next dimension of the spline"""
        spltyp = {False: -1.0, True: 0.0}
        self._features.append(feature)
        self._parameters[self.ndim*6:self.ndim*6] = (spltyp[open], low, high,
                                                     delta, lowderiv,
                                                     highderiv)
        self.ndim = self.ndim + 1

    def _parse_values(self, values, modal, param, idim):
        if len(modal) <= idim:
            modal.append(len(values))
        else:
            if modal[idim] != len(values):
                raise ValueError("Spline parameters should be a " +
                                 "rectangular matrix")
        for elmnt in values:
            if modutil.non_string_iterable(elmnt):
                self._parse_values(elmnt, modal, param, idim + 1)
            else:
                if idim + 1 != len(modal):
                    raise ValueError("Spline parameters should be a " +
                                     "rectangular matrix")
                param.append(elmnt)


# Modeller 9 compatibility
class restraint_form(RestraintForm):
    _builtin_index = RestraintForm._builtin_index

    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(restraint_form)
        RestraintForm.__init__(self, *args, **keys)


class lower_bound(LowerBound):
    _builtin_index = LowerBound._builtin_index

    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(lower_bound)
        LowerBound.__init__(self, *args, **keys)


class upper_bound(UpperBound):
    _builtin_index = UpperBound._builtin_index

    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(upper_bound)
        UpperBound.__init__(self, *args, **keys)


class gaussian(Gaussian):
    _builtin_index = Gaussian._builtin_index

    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(gaussian)
        Gaussian.__init__(self, *args, **keys)


class multi_gaussian(MultiGaussian):
    _builtin_index = MultiGaussian._builtin_index

    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(multi_gaussian)
        MultiGaussian.__init__(self, *args, **keys)


class lennard_jones(LennardJones):
    _builtin_index = LennardJones._builtin_index

    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(lennard_jones)
        LennardJones.__init__(self, *args, **keys)


class coulomb(Coulomb):
    _builtin_index = Coulomb._builtin_index

    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(coulomb)
        Coulomb.__init__(self, *args, **keys)


class cosine(Cosine):
    _builtin_index = Cosine._builtin_index

    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(cosine)
        Cosine.__init__(self, *args, **keys)


class factor(Factor):
    _builtin_index = Factor._builtin_index

    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(factor)
        Factor.__init__(self, *args, **keys)


class multi_binormal(MultiBinormal):
    _builtin_index = MultiBinormal._builtin_index

    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(multi_binormal)
        MultiBinormal.__init__(self, *args, **keys)


class spline(Spline):
    _builtin_index = Spline._builtin_index

    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(spline)
        Spline.__init__(self, *args, **keys)


class nd_spline(NDSpline):
    _builtin_index = NDSpline._builtin_index

    def __init__(self, *args, **keys):
        _deprecation_handler._init_class(nd_spline)
        NDSpline.__init__(self, *args, **keys)

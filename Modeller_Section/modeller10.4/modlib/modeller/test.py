import unittest
import modeller
import math
from modeller.util.deprecation import _deprecation_handler
try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO
import sys

__all__ = ['ModellerTest', 'TemporaryReplace']


# If we're using Python 2.6, 3.0, or 3.1, add in more modern unittest
# convenience methods
if not hasattr(unittest.TestCase, 'assertIsInstance'):
    def assertIn(self, member, container, msg=None):
        return self.assertTrue(
            member in container,
            msg or '%s not found in %s' % (member, container))

    def assertNotIn(self, member, container, msg=None):
        return self.assertTrue(
            member not in container,
            msg or '%s unexpectedly found in %s' % (member, container))

    def assertIs(self, a, b, msg=None):
        return self.assertTrue(a is b, msg or '%s is not %s' % (a, b))

    def assertIsInstance(self, obj, cls, msg=None):
        return self.assertTrue(
            isinstance(obj, cls),
            msg or '%s is not an instance of %s' % (obj, cls))

    def assertLess(self, a, b, msg=None):
        return self.assertTrue(a < b, msg or '%s not less than %s' % (a, b))

    def assertGreater(self, a, b, msg=None):
        return self.assertTrue(a > b, msg or '%s not greater than %s' % (a, b))

    def assertLessEqual(self, a, b, msg=None):
        return self.assertTrue(
            a <= b, msg or '%s not less than or equal to %s' % (a, b))

    def assertGreaterEqual(self, a, b, msg=None):
        return self.assertTrue(
            a >= b, msg or '%s not greater than or equal to %s' % (a, b))

    def assertIsNone(self, obj, msg=None):
        return self.assertTrue(obj is None, msg or '%s is not None' % obj)

    def assertIsNotNone(self, obj, msg=None):
        return self.assertTrue(obj is not None, msg or 'unexpectedly None')

    def assertAlmostEqual(self, first, second, places=None, msg=None,
                          delta=None):
        if first == second:
            return
        if delta is not None and places is not None:
            raise TypeError("specify delta or places not both")
        diff = abs(first - second)
        if delta is not None:
            if diff <= delta:
                return
            standard_msg = ("%s != %s within %s delta (%s difference)"
                            % (first, second, delta, diff))
        else:
            if places is None:
                places = 7
            if round(diff, places) == 0:
                return
            standard_msg = ("%s != %s within %r places (%s difference)"
                            % (first, second, places, diff))
        raise self.failureException(msg or standard_msg)
    unittest.TestCase.assertIn = assertIn
    unittest.TestCase.assertNotIn = assertNotIn
    unittest.TestCase.assertIs = assertIs
    unittest.TestCase.assertIsInstance = assertIsInstance
    unittest.TestCase.assertLess = assertLess
    unittest.TestCase.assertGreater = assertGreater
    unittest.TestCase.assertLessEqual = assertLessEqual
    unittest.TestCase.assertGreaterEqual = assertGreaterEqual
    unittest.TestCase.assertIsNone = assertIsNone
    unittest.TestCase.assertIsNotNone = assertIsNotNone
    unittest.TestCase.assertAlmostEqual = assertAlmostEqual


class TemporaryReplace(object):
    """Temporarily replace the named attributes. When this object is deleted,
       the attributes are restored to the original values.
       For example, to redirect stdout and stderr to a StringIO object
       for the lifetime of m::

         sio = StringIO()
         m = TemporaryReplace([(sys, 'stdout', sio), (sys, 'stderr', sio)])
    """

    def __init__(self, replacements):
        self._replacements = replacements
        self._orig_val = []
        for obj, attr, val in self._replacements:
            self._orig_val.append(getattr(obj, attr))
            setattr(obj, attr, val)

    def __del__(self):
        if hasattr(self, '_replacements'):
            for (repl, orig) in zip(self._replacements, self._orig_val):
                obj, attr, val = repl
                setattr(obj, attr, orig)


class ModellerTest(unittest.TestCase):
    _env = None

    def get_environ(self):
        """Set up the Modeller environment, and keep a reference to it so that
           we only need to do it once for all tests"""
        if not self._env:
            env = modeller.Environ()
            env.io.atom_files_directory = ['../data', 'test/data']
            env.edat.dynamic_sphere = False
            env.libs.topology.read(file='${LIB}/top_heav.lib')
            env.libs.parameters.read(file='${LIB}/par.lib')
            ModellerTest._env = env
        return self._env

    def assertModelsEqual(self, mdl1, mdl2):
        """Check to make sure that models ``mdl1`` and ``mdl2`` are the same"""
        self.assertEqual(len(mdl1.atoms), len(mdl2.atoms))
        for (at1, at2) in zip(mdl1.atoms, mdl2.atoms):
            self.assertAtomsEqual(at1, at2)

    def assertAtomsEqual(self, at1, at2):
        """Check to make sure atoms ``at1`` and ``at2`` are the same"""
        self.assertEqual(at1.name, at2.name)
        self.assertAlmostEqual(at1.x, at2.x, places=3)

    def assertDistance(self, o1, o2, mean, tolerance):
        """Check to make sure the distance between o1 and o2 is as expected"""
        distsq = (o1.x-o2.x)**2 + (o1.y-o2.y)**2 + (o1.z-o2.z)**2
        dist = distsq ** 0.5
        msg = "Distance between %s and %s is %f - expected %f" % (o1, o2,
                                                                  dist, mean)
        self.assertAlmostEqual(dist, mean, msg=msg, delta=tolerance)

    def assertAlignmentsEqual(self, refaln, aln, check_meta=True):
        """Check to make sure that alignments ``refaln`` and ``aln``
           are the same"""
        self.assertEqual(len(refaln), len(aln),
                         "Inconsistent number of sequences (%d vs. %d)"
                         % (len(refaln), len(aln)))
        for (seq1, seq2) in zip(refaln, aln):
            self.assertAlnSequencesEqual(seq1, seq2, check_meta)
        self.assertEqual(len(refaln.comments), len(aln.comments))
        for (c1, c2) in zip(refaln.comments, aln.comments):
            self.assertEqual(c1, c2)

    def assertAlnSequencesEqual(self, seq1, seq2, check_meta=True):
        """Check to make sure that sequences ``seq1`` and ``seq2``
           are the same"""
        self.assertEqual(len(seq1), len(seq2),
                         "Inconsistent number of residues (%d vs. %d)"
                         % (len(seq1), len(seq2)))
        if check_meta:
            self.assertEqual(seq1.range[0], seq2.range[0])
            self.assertEqual(seq1.range[1], seq2.range[1])
            self.assertEqual(seq1.resolution, seq2.resolution)
            self.assertEqual(seq1.rfactor, seq2.rfactor)
            self.assertEqual(seq1.code, seq2.code)
            self.assertEqual(seq1.atom_file, seq2.atom_file)
            self.assertEqual(seq1.name, seq2.name)
            self.assertEqual(seq1.prottyp, seq2.prottyp)
            self.assertEqual(seq1.source, seq2.source)
        for (res1, res2) in zip(seq1.residues, seq2.residues):
            self.assertResiduesEqual(res1, res2)
            self.assertEqual(res1.get_position().num, res2.get_position().num)
            self.assertEqual(res1.get_leading_gaps(), res2.get_leading_gaps())
            self.assertEqual(res1.get_trailing_gaps(),
                             res2.get_trailing_gaps())

    def assertResiduesEqual(self, res1, res2):
        """Check to make sure that residues ``res1`` and ``res2``
           are the same"""
        self.assertEqual(res1.name, res2.name,
                         "Inconsistent residue names (%s vs. %s)"
                         % (res1.name, res2.name))

    def assertProfilesEqual(self, prf1, prf2):
        """Check to make sure that two profile objects are the same"""
        self.assertEqual(len(prf1), len(prf2))
        for (seq1, seq2) in zip(prf1, prf2):
            self.assertProfileSeqsEqual(seq1, seq2)

    def assertProfileSeqsEqual(self, seq1, seq2):
        """Check to make sure that two profile sequences are the same"""
        self.assertEqual(len(seq1), len(seq2))
        self.assertEqual(seq1.code, seq2.code)
        self.assertEqual(seq1.prottyp, seq2.prottyp)
        self.assertEqual(seq1.iter, seq2.iter)
        self.assertEqual(seq1.neqv, seq2.neqv)
        self.assertEqual(seq1.fid, seq2.fid)
        self.assertEqual(seq1.evalue, seq2.evalue)

    def assertSeqDBsEqual(self, sdb1, sdb2):
        """Check to make sure that two sequence_db objects are the same"""
        self.assertEqual(len(sdb1), len(sdb2))
        for (seq1, seq2) in zip(sdb1, sdb2):
            self.assertSeqDBseqsEqual(seq1, seq2)

    def assertSeqDBseqsEqual(self, seq1, seq2):
        """Check to make sure that two sequence_db sequences are the same"""
        self.assertEqual(len(seq1), len(seq2))
        self.assertEqual(seq1.code, seq2.code)
        self.assertEqual(seq1.prottyp, seq2.prottyp)
        self.assertEqual(seq1.resol, seq2.resol)
        self.assertEqual(len(seq1.residues), len(seq2.residues))
        for (r1, r2) in zip(seq1.residues, seq2.residues):
            self.assertSeqDBresEqual(r1, r2)

    def assertSeqDBresEqual(self, res1, res2):
        """Check to make sure that two sequence_db residues are the same"""
        self.assertEqual(res1.type, res2.type)

    def run_capture_stdout(self, method, *args, **keys):
        """Run a method and capture its standard output. Returns both the
           method's own return values and the standard output."""
        saved_stdout = sys.stdout
        sio = StringIO()
        try:
            sys.stdout = sio
            ret = method(*args, **keys)
        finally:
            sys.stdout = saved_stdout
        val = sio.getvalue()
        sys.stdout.write(val)
        return ret, val

    def assertIsDeprecatedAlias(self, depcls, cls, *args, **keys):
        """Make sure that class `depcls` is a deprecated alias for `cls`"""
        create_alias = keys.pop('create_alias', True)
        self.assertTrue(issubclass(depcls, cls))
        old_depmode = _deprecation_handler.depmode
        try:
            if create_alias:
                # try creating the deprecated class; make sure we get a warning
                _deprecation_handler.depmode = "WARN"
                ret, output = self.run_capture_stdout(depcls, *args, **keys)
                self.assertIn('is deprecated; use', output)
            else:
                ret = None
            # try creating the deprecated class; make sure we get an error
            _deprecation_handler.depmode = "ERROR"
            self.assertRaises(modeller.ModellerError, depcls, *args, **keys)
            return ret
        finally:
            _deprecation_handler.depmode = old_depmode

    def create_angle(self, angle, a1, a2, a3):
        """Set the coordinates of the given atoms to form (by construction)
           the given angle between them."""
        a1.x, a1.y, a1.z = -10.0, 0.0, 0.0
        a2.x, a2.y, a2.z = 0.0, 0.0, 0.0
        a3.x, a3.y, a3.z = -math.cos(angle), math.sin(angle), 0.0

    def create_dihedral(self, angle, a1, a2, a3, a4):
        """Set the coordinates of the given atoms to form (by construction)
           the given dihedral angle between them."""
        a1.x, a1.y, a1.z = 10.0, 0.0, -10.0
        a2.x, a2.y, a2.z = 0.0, 0.0, -10.0
        a3.x, a3.y, a3.z = 0.0, 0.0, 0.0
        a4.x, a4.y, a4.z = math.cos(angle), math.sin(angle), 0.0

    def raw_byte_string(self, s):
        """Return a raw byte (not Unicode) string constant"""
        if sys.version_info[0] >= 3:
            # Cannot use b'foo' since that is invalid syntax in Python 2.3
            return bytes(s, encoding='ascii')
        else:
            return s

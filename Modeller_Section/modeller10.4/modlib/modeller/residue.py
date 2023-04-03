"""Classes for handling residues."""

import _modeller

__docformat__ = "epytext en"


class Residue(object):
    """A single residue in a sequence or structure"""

    def __init__(self, mdl, num):
        self.mdl = mdl
        self._num = num

    def __eq__(self, r):
        # Two Residue objects are considered equal iff they represent the
        # same residue in the same sequence
        return type(r) is type(self) and self.mdl is r.mdl \
               and self._num == r._num

    def __ne__(self, r):
        return not self == r

    def __hash__(self):
        return hash((self.mdl, self._num))

    def __repr__(self):
        return "Residue %s" % self.name

    def __str__(self):
        return "<%s>" % repr(self)

    def __get_name(self):
        return _modeller.mod_residue_name_from_type(self.type,
                                                    self.mdl.env.libs.modpt)

    def __set_name(self, val):
        t = _modeller.mod_residue_type_from_name(val, self.mdl.env.libs.modpt)
        if t == 0:
            raise ValueError("Invalid residue name: '%s'" % str(val))
        self.type = t

    def __get_pdb_name(self):
        return _modeller.mod_residue_pdbnam_from_type(self.type,
                                                      self.mdl.env.libs.modpt)

    def __set_pdb_name(self, val):
        t = _modeller.mod_residue_type_from_pdbnam(val,
                                                   self.mdl.env.libs.modpt)
        if t == 0:
            raise ValueError("Invalid PDB residue name: '%s'" % str(val))
        self.type = t

    def __get_code(self):
        return _modeller.mod_residue_code_from_type(self.type,
                                                    self.mdl.env.libs.modpt)

    def __set_code(self, val):
        t = _modeller.mod_residue_type_from_code(val, self.mdl.env.libs.modpt)
        if t == 0:
            raise ValueError("Invalid residue code: '%s'" % str(val))
        self.type = t

    def __get_hetatm(self):
        return _modeller.mod_residue_is_hetatm(self.type,
                                               self.mdl.env.libs.modpt)

    name = property(__get_name, __set_name,
                    doc="Internal (CHARMM) residue type name")
    pdb_name = property(__get_pdb_name, __set_pdb_name,
                        doc="PDB/mmCIF (IUPAC) type name")
    code = property(__get_code, __set_code, doc="One-letter residue type code")
    hetatm = property(__get_hetatm, doc="Whether this is a PDB HETATM residue")

from modeller.alignment import Alignment
from modeller.selection import Selection
from modeller.model import Model


def fit(env, model, code, model2, code2, alnfile, model2_fit):
    """Superposes model2 on model, and writes out a file with model2 superposed
       on model."""
    m1 = Model(env, file=model)
    m2 = Model(env, file=model2)

    aln = Alignment(env, file=alnfile, align_codes=(code, code2))
    atmsel = Selection(m1).only_atom_types('CA')

    atmsel.superpose(m2, aln)

    m2.write(file=model2_fit)

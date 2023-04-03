from modeller import Model, Alignment


def complete_pdb(env, filename, special_patches=None, transfer_res_num=False,
                 model_segment=None, patch_default=True):
    """Reads the given PDB file, reorders the atoms to match the current
       topology library, and adds any missing atoms.

       You should read topology and parameters into 'env' before calling
       this routine.

       :param env: Modeller environment.
       :type  env: :class:`Environ`
       :param filename: the PDB file to read.
       :param special_patches: if set, it is expected to be a routine which
              takes one parameter (the model) and applies any patches (e.g.
              disulfide bridges).
       :param transfer_res_num: if True, the residue numbering from the
              original PDB is retained (by default, residues are renumbered
              from 1).
       :param patch_default: if True, default terminal patches are applied.

       :return: the completed model.
       :rtype: :class:`Model`"""

    vars = {}
    if model_segment is not None:
        vars['model_segment'] = model_segment
    mdl = Model(env, file=filename, model_format='PDB_OR_MMCIF', **vars)
    # Save original chain IDs, since generate_topology resets them
    chain_ids = [c.name for c in mdl.chains]
    aln = Alignment(env)
    aln.append_model(mdl, atom_files=filename, align_codes='struc')
    aln.append_model(mdl, atom_files=filename+'.ini', align_codes='struc-ini')
    mdl.clear_topology()
    mdl.generate_topology(aln[-1], patch_default=patch_default)
    if special_patches:
        special_patches(mdl)
    # Save original seq_id, as transfer_xyz sets it
    seq_id = mdl.seq_id
    mdl.transfer_xyz(aln)
    mdl.seq_id = seq_id
    # Restore original chain IDs
    for (chain, chainid) in zip(mdl.chains, chain_ids):
        chain.name = chainid
    mdl.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')
    if transfer_res_num:
        mdl2 = Model(env, file=filename, **vars)
        mdl.res_num_from(mdl2, aln)
    return mdl

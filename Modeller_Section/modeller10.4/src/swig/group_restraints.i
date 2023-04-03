struct mod_atom_classes {
%immutable;
  int ngrpatm;
  int nattmod;
  struct mod_int1_array iattmod;
};

struct mod_group_restraints {
%immutable;
  struct mod_atom_classes atclass;
};

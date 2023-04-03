%apply char **STRVEC { char **atom_files_directory };
%typemap(freearg) (char **atom_files_directory) {
  /* Don't free atom_files_directory, since mod_io_data_set uses it */
}

%inline %{

/** Set members of an io_data object. */
static void mod_io_data_set(struct mod_io_data *io, gboolean hydrogen,
                            gboolean hetatm, gboolean water,
                            gboolean convert_modres,
                            gboolean hybrid36, gboolean two_char_chain,
                            char **atom_files_directory)
{
  io->hydrogen = hydrogen;
  io->hetatm = hetatm;
  io->water = water;
  io->convert_modres = convert_modres;
  io->hybrid36 = hybrid36;
  io->two_char_chain = two_char_chain;
  g_strfreev(io->atom_files_directory);
  /* Note that we don't copy the vector here,
     so SWIG must NOT free the input! */
  io->atom_files_directory = atom_files_directory;
}

%}

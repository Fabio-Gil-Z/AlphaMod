struct mod_coordinates {
%immutable;
  int natm;
  int dirty;
  struct mod_int1_array iresatm;
  struct mod_int1_array iatmr1;
  struct mod_float1_array occ;
  struct mod_float1_array biso;
  struct mod_float1_array atmacc;
  struct mod_float1_array x, y, z;
};

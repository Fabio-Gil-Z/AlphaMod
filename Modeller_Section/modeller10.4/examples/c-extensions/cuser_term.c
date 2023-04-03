#include <glib.h>
#include <stdlib.h>
#include "modeller.h"
#include "cuser_term.h"

/** Data used by custom energy term */
struct myterm_data {
  float strength;
};

/** Evaluate the custom energy term */
static int myterm_eval(void *data, const struct mod_model *model,
                       gboolean deriv, const int *atind, int n_atind,
                       float *e_term, float *dvx, float *dvy, float *dvz,
                       struct mod_profile_data *profile,
                       const struct mod_energy_data *enedata,
                       const struct mod_libraries *libs)
{
  int i;
  struct myterm_data *mytd = (struct myterm_data *)data;
  *e_term = 0.0;
  if (deriv) {
    for (i = 0; i < n_atind; i++) {
      dvx[i] = dvy[i] = dvz[i] = 0.0;
    }
  }
  for (i = 0; i < n_atind; i++) {
    int iatm = atind[i] - 1;
    float x = mod_float1_get(&model->cd.x, iatm);
    if (x > 10.0) {
      float e = (x - 10.0) * mytd->strength;
      *e_term += e;
      if (profile) {
        mod_energy_term_add_profile(profile, &iatm, 1, e, 0., 0., 0.);
      }
      if (deriv) {
        dvx[i] += mytd->strength;
      }
    }
  }
  return 0;
}

/** Create the custom energy term */
void myterm_create(struct mod_energy_data *edat, int indx, int physical_type,
                   float strength)
{
  struct myterm_data *mytd = malloc(sizeof(struct myterm_data));
  mytd->strength = strength;
  mod_energy_term_new(edat, indx, myterm_eval, free, mytd, 0.0, physical_type);
}

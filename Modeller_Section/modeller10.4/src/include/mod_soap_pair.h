/** \file mod_soap_pair.h           Functions for SOAP pairwise scoring.
 *
 *             Part of MODELLER, Copyright(c) 1989-2022 Andrej Sali
 */

#ifndef MOD_SOAP_PAIR_H
#define MOD_SOAP_PAIR_H

#include "mod_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Create the SOAP pairwise energy term */
void mod_soap_pair_create(struct mod_energy_data *edat, int indx,
                          int physical_type, const char *library,
                          int *ierr);

#ifdef __cplusplus
}
#endif
#endif  /* MOD_SOAP_PAIR_H */

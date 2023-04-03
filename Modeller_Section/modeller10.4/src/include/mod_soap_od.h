/** \file mod_soap_od.h        Functions for orientation-dependent SOAP scoring.
 *
 *             Part of MODELLER, Copyright(c) 1989-2022 Andrej Sali
 */

#ifndef MOD_SOAP_OD_H
#define MOD_SOAP_OD_H

#include "mod_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Create the orientation-dependent SOAP energy term */
void mod_soap_od_create(struct mod_energy_data *edat, int indx,
                        int physical_type, const char *library, int *ierr);

#ifdef __cplusplus
}
#endif
#endif  /* MOD_SOAP_OD_H */

/** \file mod_soap_access.h           Functions for SOAP accessiblity scoring.
 *
 *             Part of MODELLER, Copyright(c) 1989-2022 Andrej Sali
 */

#ifndef MOD_SOAP_ACCESS_H
#define MOD_SOAP_ACCESS_H

#include "mod_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Create the SOAP atom accessiblity energy term */
void mod_soap_access_create(struct mod_energy_data *edat, int indx,
                            int physical_type, const char *library);

#ifdef __cplusplus
}
#endif
#endif  /* MOD_SOAP_ACCESS_H */

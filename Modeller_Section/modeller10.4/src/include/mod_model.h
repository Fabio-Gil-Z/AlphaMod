/** \file mod_model.h        Functions for handling mod_model objects.
 *
 *             Part of MODELLER, Copyright(c) 1989-2022 Andrej Sali
 */

#ifndef MOD_MODEL_H
#define MOD_MODEL_H

#include <glib.h>
#include "mod_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Store some data in the model, indexed by the given key.
    When the data is released, freefunc is called (if is is non-NULL) to
    release any dynamic memory associated with it.
    The key is usually a pointer to another data structure; note that
    mod_model_release_data() should be called when that data structure is
    freed.
 */
void mod_model_set_data(const struct mod_model *model, gpointer key,
                        gpointer data, GDestroyNotify freefunc);

/** Retrieve data from the model indexed by the key.
    If no data is available, NULL is returned (note that it is not possible
    to distinguish this from the data being 'NULL').
 */
gpointer mod_model_get_data(const struct mod_model *model, gpointer key);

/** Release any data indexed by the key from all models. */
void mod_model_release_data(gpointer key);

#ifdef __cplusplus
}
#endif
#endif  /* MOD_MODEL_H */

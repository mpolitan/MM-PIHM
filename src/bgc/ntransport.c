#include "pihm.h"

void NTransport (pihm_struct pihm)
{
    int             i;

    /*
     * Calculate solute N concentrantions
     */
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        elem_struct *elem;
        int         j;
        double      strg;

        elem = &pihm->elem[i];

        /* Initialize N fluxes */
        for (j = 0; j < NUM_EDGE; j++)
        {
            elem->nsol.subflux[j] = 0.0;
        }

        /* Element subsurface */
        strg = (elem->ws.unsat + elem->ws.gw) * elem->soil.porosity +
            elem->soil.depth * elem->soil.smcmin;
        elem->nsol.conc_subsurf = (strg > 0.0) ?
            elem->ns.sminn / strg / 1000.0 : 0.0;
        elem->nsol.conc_subsurf = (elem->nsol.conc_subsurf > 0.0) ?
            elem->nsol.conc_subsurf : 0.0;
    }

    /*
     * Calculate solute fluxes
     */
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        elem_struct *elem;
        int         j;

        elem = &pihm->elem[i];

        /* Element to element */
        for (j = 0; j < NUM_EDGE; j++)
        {
            if (elem->nabr[j] > 0)
            {
                elem->nsol.subflux[j] = elem->wf.subsurf[j] * 1000.0 *
                    ((elem->wf.subsurf[j] > 0.0) ?
                    MOBILEN_PROPORTION * elem->nsol.conc_subsurf : 0.0);
            }
            else
            {
                elem->nsol.subflux[j] = 0.0;
            }
        }
    }
}

#include "pihm.h"

void LateralFlow (pihm_struct pihm)
{
    int             i;
    double         *dhbydx;
    double         *dhbydy;

    dhbydx = (double *)malloc (nelem * sizeof (double));
    dhbydy = (double *)malloc (nelem * sizeof (double));

    FrictSlope (pihm->elem, pihm->riv, pihm->ctrl.surf_mode, dhbydx, dhbydy);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
        int         j;
        double      dif_y_sub;
        double      avg_y_sub;
        double      grad_y_sub;
        double      effk;
        double      effk_nabr;
        double      avg_ksat;
        double      dif_y_surf;
        double      avg_y_surf;
        double      grad_y_surf;
        double      avg_sf;
        double      avg_rough;
        double      crossa;

        elem_struct *elem;
        elem_struct *nabr;

        elem = &pihm->elem[i];

        for (j = 0; j < NUM_EDGE; j++)
        {
            if (elem->nabr[j] > 0)
            {
                nabr = &pihm->elem[elem->nabr[j] - 1];

                /*
                 * Subsurface lateral flux calculation between triangular
                 * elements
                 */
                dif_y_sub =
                    (elem->ws.gw + elem->topo.zmin) - (nabr->ws.gw +
                    nabr->topo.zmin);
                avg_y_sub = AvgY (dif_y_sub, elem->ws.gw, nabr->ws.gw);
                grad_y_sub = dif_y_sub / elem->topo.nabrdist[j];
                /* Take into account macropore effect */
                effk =
                    EffKH (elem->ws.gw, elem->soil.depth, elem->soil.dmac,
                    elem->soil.kmach, elem->soil.areafv, elem->soil.ksath);
                effk_nabr =
                    EffKH (nabr->ws.gw, nabr->soil.depth, nabr->soil.dmac,
                    nabr->soil.kmach, nabr->soil.areafv, nabr->soil.ksath);
                avg_ksat = 0.5 * (effk + effk_nabr);
                /* Groundwater flow modeled by Darcy's Law */
                elem->wf.subsurf[j] =
                    avg_ksat * grad_y_sub * avg_y_sub * elem->topo.edge[j];

                /*
                 * Surface lateral flux calculation between triangular
                 * elements
                 */
                if (pihm->ctrl.surf_mode == KINEMATIC)
                {
                    dif_y_surf = elem->topo.zmax - nabr->topo.zmax;
                }
                else
                {
                    dif_y_surf = (elem->ws.surfh + elem->topo.zmax) -
                        (nabr->ws.surfh + nabr->topo.zmax);
                }
                avg_y_surf = AvgYsfc (dif_y_surf, elem->ws.surfh,
                    nabr->ws.surfh);
                grad_y_surf = dif_y_surf / elem->topo.nabrdist[j];
                avg_sf = 0.5 *
                    (sqrt (dhbydx[i] * dhbydx[i] + dhbydy[i] * dhbydy[i]) +
                    sqrt (dhbydx[nabr->ind - 1] * dhbydx[nabr->ind - 1] +
                    dhbydy[nabr->ind - 1] * dhbydy[nabr->ind - 1]));
                if (pihm->ctrl.surf_mode == KINEMATIC)
                {
                    avg_sf = (grad_y_surf > 0.0) ? grad_y_surf : GRADMIN;
                }
                else
                {
                    avg_sf = (avg_sf > GRADMIN) ? avg_sf : GRADMIN;
                }
                /* Weighting needed */
                avg_rough = 0.5 * (elem->lc.rough + nabr->lc.rough);
                crossa = avg_y_surf * elem->topo.edge[j];
                elem->wf.ovlflow[j] =
                    OverlandFlow (avg_y_surf, grad_y_surf, avg_sf, crossa,
                    avg_rough);
            }
            else if (elem->nabr[j] < 0)
            {
                /* Do nothing. River-element interactions are calculated
                 * in river_flow.c */
            }
            else                        /* Boundary condition flux */
            {
                /* No flow (natural) boundary condition is default */
                if (elem->attrib.bc_type[j] == 0)
                {
                    elem->wf.ovlflow[j] = 0.0;
                    elem->wf.subsurf[j] = 0.0;
                }
                /* Note: ideally different boundary conditions need to be
                 * incorporated for surf and subsurf respectively */
                else if (elem->attrib.bc_type[j] > 0)
                {
                    /* Note: the formulation assumes only Dirichlet TS right
                     * now */
                    /* note the assumption here is no flow for surface */
                    elem->wf.ovlflow[j] = 0.0;
                    dif_y_sub =
                        elem->ws.gw + elem->topo.zmin - elem->bc.head[j];
                    avg_y_sub =
                        AvgY (dif_y_sub, elem->ws.gw,
                        elem->bc.head[j] - elem->topo.zmin);
                    /* Minimum distance from circumcenter to the edge of the
                     * triangle on which boundary condition is defined */
                    effk =
                        EffKH (elem->ws.gw, elem->soil.depth, elem->soil.dmac,
                        elem->soil.kmach, elem->soil.areafv,
                        elem->soil.ksath);
                    avg_ksat = effk;
                    grad_y_sub = dif_y_sub / elem->topo.nabrdist[j];
                    elem->wf.subsurf[j] =
                        avg_ksat * grad_y_sub * avg_y_sub *
                        elem->topo.edge[j];
                }
                else
                {
                    /* Neumann bc (note: md->ele[i].bc[j] value has to be
                     * = 2+(index of neumann boundary ts) */
                    elem->wf.ovlflow[j] = elem->bc.flux[j];
                    elem->wf.subsurf[j] = elem->bc.flux[j];
                }
            }                   /* End of specified boundary condition */
        }                       /* End of neighbor loop */
    }                           /* End of element loop */

    free (dhbydx);
    free (dhbydy);

#ifdef _FBR_
    FBRLatFlow(pihm->elem, pihm->riv);
#endif
}

void FrictSlope (elem_struct *elem, river_struct *riv, int surf_mode,
    double *dhbydx, double *dhbydy)
{
    int             i;
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (i = 0; i < nelem; i++)
    {
    int             j;
    double          surfh[NUM_EDGE];
    elem_struct    *nabr;
    river_struct   *rivnabr;

        if (surf_mode == DIFF_WAVE)
        {
            for (j = 0; j < NUM_EDGE; j++)
            {
                if (elem[i].nabr[j] > 0)
                {
                    nabr = &elem[elem[i].nabr[j] - 1];
                    surfh[j] = nabr->topo.zmax + nabr->ws.surfh;
                }
                else if (elem[i].nabr[j] < 0)
                {
                    rivnabr = &riv[-elem[i].nabr[j] - 1];

                    if (rivnabr->ws.stage > rivnabr->shp.depth)
                    {
                        surfh[j] = rivnabr->topo.zbed + rivnabr->ws.stage;
                    }
                    else
                    {
                        surfh[j] = rivnabr->topo.zmax;
                    }
                }
                else
                {
                    if (elem[i].attrib.bc_type[j] == 0)
                    {
                        surfh[j] = elem[i].topo.zmax + elem[i].ws.surfh;
                    }
                    else
                    {
                        surfh[j] = elem[i].bc.head[j];
                    }
                }
            }

            dhbydx[i] = DhByDl (elem[i].topo.nabrdist_y,
                elem[i].topo.nabrdist_x, surfh);
            dhbydy[i] = DhByDl (elem[i].topo.nabrdist_x,
                elem[i].topo.nabrdist_y, surfh);
        }
    }
}

double AvgYsfc (double diff, double yi, double yinabr)
{
    double          avg_y;

    if (diff > 0.0)
    {
        if (yi > DEPRSTG)
        {
            avg_y = 1.0 * (yi - DEPRSTG);
        }
        else
        {
            avg_y = 0.0;
        }
    }
    else
    {
        if (yinabr > DEPRSTG)
        {
            avg_y = 1.0 * (yinabr - DEPRSTG);
        }
        else
        {
            avg_y = 0.0;
        }
    }

    return (avg_y);
}

double AvgY (double diff, double yi, double yinabr)
{
    double          avg_y;

    if (diff > 0.0)
    {
        if (yi > 0.0)
        {
            avg_y = 1.0 * yi;
        }
        else
        {
            avg_y = 0.0;
        }
    }
    else
    {
        if (yinabr > 0.0)
        {
            avg_y = 1.0 * yinabr;
        }
        else
        {
            avg_y = 0.0;
        }
    }

    return (avg_y);
}

double DhByDl (double *l1, double *l2, double *surfh)
{
    return (-1.0 * (l1[2] * (surfh[1] - surfh[0]) +
        l1[1] * (surfh[0] - surfh[2]) + l1[0] * (surfh[2] - surfh[1])) /
        (l2[2] * (l1[1] - l1[0]) + l2[1] * (l1[0] - l1[2]) +
        l2[0] * (l1[2] - l1[1])));
}

double EffKH (double tmpy, double aqdepth, double macd, double macksath,
    double areaf, double ksath)
{
    tmpy = (tmpy > 0.0) ? tmpy : 0.0;

    if (tmpy > aqdepth - macd)
    {
        if (tmpy > aqdepth)
        {
            return (macksath * macd * areaf + ksath * (aqdepth -
                    macd * areaf)) / aqdepth;
        }
        else
        {
            return (macksath * (tmpy - (aqdepth - macd)) * areaf +
                ksath * (aqdepth - macd + (tmpy - (aqdepth -
                            macd)) * (1.0 - areaf))) / tmpy;
        }
    }
    else
    {
        return (ksath);
    }
}


double OverlandFlow (double avg_y, double grad_y, double avg_sf,
    double crossa, double avg_rough)
{
    return (crossa * pow (avg_y,
            2.0 / 3.0) * grad_y / (sqrt (fabs (avg_sf)) * avg_rough));
}

#ifdef _FBR_
void FBRLatFlow (elem_struct *elem, river_struct *riv)
{
    /*
     * Fractured bedrock lateral flux calculation between triangular
     * elements
     */
    int             i, j;
    elem_struct    *nabr;
    double          dif_y_fbr;
    double          avg_y_fbr;
    double          grad_y_fbr;
    double          avg_ksat;

    for (i = 0; i < nelem; i++)
    {
        for (j = 0; j < NUM_EDGE; j++)
        {
            if (elem[i].nabr[j] == 0)
            {
                elem[i].wf.fbrflow[j] = 0.0;
            }
            else
            {
                if (elem[i].nabr[j] > 0)
                {
                    nabr = &elem[elem[i].nabr[j]];
                }
                else
                {
                    nabr = (riv[-elem[i].nabr[j] - 1].leftele == elem[i].ind) ?
                        &elem[riv[-elem[i].nabr[j] - 1].rightele - 1] :
                        &elem[riv[-elem[i].nabr[j] - 1].leftele - 1];
                }

                dif_y_fbr =
                    (elem->ws.fbr_gw + elem->topo.zbed) -
                    (nabr->ws.fbr_gw + nabr->topo.zbed);
                avg_y_fbr = AvgY (dif_y_fbr, elem->ws.fbr_gw, nabr->ws.fbr_gw);
                grad_y_fbr = dif_y_fbr / elem->topo.nabrdist[j];
                avg_ksat = 0.5 * (elem->geol.ksath + nabr->geol.ksath);
                /* Groundwater flow modeled by Darcy's Law */
                elem->wf.fbrflow[j] =
                    avg_ksat * grad_y_fbr * avg_y_fbr * elem->topo.edge[j];
            }
        }
    }
}
#endif

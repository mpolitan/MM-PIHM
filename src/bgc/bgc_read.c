#include "pihm.h"

void ReadBgc (char *fn, ctrl_struct *ctrl, co2control_struct *co2,
    ndepcontrol_struct *ndepctrl, cninit_struct *cninit, char *co2_fn,
    char *ndep_fn)
{
    FILE           *bgc_file;
    pihm_t_struct   pihm_time1, pihm_time2;
    char            cmdstr[MAXSTRING];
    int             lno = 0;

    /* Read bgc simulation control file */
    bgc_file = fopen (fn, "r");
    CheckFile (bgc_file, fn);
    PIHMprintf (VL_VERBOSE, " Reading %s\n", fn);

    FindLine (bgc_file, "RESTART", &lno, fn);
    NextLine (bgc_file, cmdstr, &lno);
    sscanf (cmdstr, "%d", &ctrl->read_bgc_restart);
    NextLine (bgc_file, cmdstr, &lno);
    sscanf (cmdstr, "%d", &ctrl->write_bgc_restart);

    FindLine (bgc_file, "TIME_DEFINE", &lno, fn);
    NextLine (bgc_file, cmdstr, &lno);
    sscanf (cmdstr, "%d", &spinup_mode);
    NextLine (bgc_file, cmdstr, &lno);
    sscanf (cmdstr, "%d", &ctrl->maxspinyears);

    /* In spinup mode, simulation time should be full years */
    if (spinup_mode)
    {
        pihm_time1 = PIHMTime (ctrl->starttime);
        pihm_time2 = PIHMTime (ctrl->endtime);

        if (pihm_time1.month != pihm_time2.month ||
            pihm_time1.day != pihm_time2.day ||
            pihm_time1.hour != pihm_time2.hour ||
            pihm_time1.minute != pihm_time2.minute)
        {
            PIHMprintf (VL_ERROR,
                "Error: In BGC spinup mode, "
                "simulation period should be full years.\n");
            PIHMprintf (VL_ERROR,
                "Please check your .para input file.\n");
            PIHMexit (EXIT_FAILURE);
        }
    }

    FindLine (bgc_file, "CO2_CONTROL", &lno, fn);
    NextLine (bgc_file, cmdstr, &lno);
    sscanf (cmdstr, "%d", &co2->varco2);
    NextLine (bgc_file, cmdstr, &lno);
    sscanf (cmdstr, "%lf", &co2->co2ppm);
    NextLine (bgc_file, cmdstr, &lno);
    sscanf (cmdstr, "%s", co2_fn);

    FindLine (bgc_file, "NDEP_CONTROL", &lno, fn);
    NextLine (bgc_file, cmdstr, &lno);
    sscanf (cmdstr, "%d", &ndepctrl->varndep);
    NextLine (bgc_file, cmdstr, &lno);
    sscanf (cmdstr, "%lf", &ndepctrl->ndep);
    NextLine (bgc_file, cmdstr, &lno);
    sscanf (cmdstr, "%lf", &ndepctrl->nfix);
    NextLine (bgc_file, cmdstr, &lno);
    sscanf (cmdstr, "%s", ndep_fn);

    FindLine (bgc_file, "C_STATE", &lno, fn);
    NextLine (bgc_file, cmdstr, &lno);
    sscanf (cmdstr, "%lf", &cninit->max_leafc);
    NextLine (bgc_file, cmdstr, &lno);
    sscanf (cmdstr, "%lf", &cninit->max_stemc);
    NextLine (bgc_file, cmdstr, &lno);
    sscanf (cmdstr, "%lf", &cninit->cwdc);
    NextLine (bgc_file, cmdstr, &lno);
    sscanf (cmdstr, "%lf", &cninit->litr1c);
    NextLine (bgc_file, cmdstr, &lno);
    sscanf (cmdstr, "%lf", &cninit->litr2c);
    NextLine (bgc_file, cmdstr, &lno);
    sscanf (cmdstr, "%lf", &cninit->litr3c);
    NextLine (bgc_file, cmdstr, &lno);
    sscanf (cmdstr, "%lf", &cninit->litr4c);
    NextLine (bgc_file, cmdstr, &lno);
    sscanf (cmdstr, "%lf", &cninit->soil1c);
    NextLine (bgc_file, cmdstr, &lno);
    sscanf (cmdstr, "%lf", &cninit->soil2c);
    NextLine (bgc_file, cmdstr, &lno);
    sscanf (cmdstr, "%lf", &cninit->soil3c);
    NextLine (bgc_file, cmdstr, &lno);
    sscanf (cmdstr, "%lf", &cninit->soil4c);

    FindLine (bgc_file, "N_STATE", &lno, fn);
    NextLine (bgc_file, cmdstr, &lno);
    sscanf (cmdstr, "%lf", &cninit->litr1n);
    NextLine (bgc_file, cmdstr, &lno);
    sscanf (cmdstr, "%lf", &cninit->sminn);

    FindLine (bgc_file, "DAILY_OUTPUT", &lno, fn);

    NextLine (bgc_file, cmdstr, &lno);
    ctrl->prtvrbl[LAI_CTRL] = ReadPrtCtrl (cmdstr, "LAI", fn, lno);

    NextLine (bgc_file, cmdstr, &lno);
    ctrl->prtvrbl[VEGC_CTRL] = ReadPrtCtrl (cmdstr, "VEGC", fn, lno);

    NextLine (bgc_file, cmdstr, &lno);
    ctrl->prtvrbl[LITRC_CTRL] = ReadPrtCtrl (cmdstr, "LITRC", fn, lno);

    NextLine (bgc_file, cmdstr, &lno);
    ctrl->prtvrbl[SOILC_CTRL] = ReadPrtCtrl (cmdstr, "SOILC", fn, lno);

    NextLine (bgc_file, cmdstr, &lno);
    ctrl->prtvrbl[TOTALC_CTRL] = ReadPrtCtrl (cmdstr, "TOTALC", fn, lno);

    NextLine (bgc_file, cmdstr, &lno);
    ctrl->prtvrbl[NPP_CTRL] = ReadPrtCtrl (cmdstr, "NPP", fn, lno);

    NextLine (bgc_file, cmdstr, &lno);
    ctrl->prtvrbl[NEP_CTRL] = ReadPrtCtrl (cmdstr, "NEP", fn, lno);

    NextLine (bgc_file, cmdstr, &lno);
    ctrl->prtvrbl[NEE_CTRL] = ReadPrtCtrl (cmdstr, "NEE", fn, lno);

    NextLine (bgc_file, cmdstr, &lno);
    ctrl->prtvrbl[GPP_CTRL] = ReadPrtCtrl (cmdstr, "GPP", fn, lno);

    NextLine (bgc_file, cmdstr, &lno);
    ctrl->prtvrbl[SMINN_CTRL] = ReadPrtCtrl (cmdstr, "SMINN", fn, lno);

    NextLine (bgc_file, cmdstr, &lno);
    ctrl->prtvrbl[LEAFC_CTRL] = ReadPrtCtrl (cmdstr, "LEAFC", fn, lno);

    NextLine (bgc_file, cmdstr, &lno);
    ctrl->prtvrbl[LIVESTEMC_CTRL] =
        ReadPrtCtrl (cmdstr, "LIVESTEMC", fn, lno);

    NextLine (bgc_file, cmdstr, &lno);
    ctrl->prtvrbl[DEADSTEMC_CTRL] =
        ReadPrtCtrl (cmdstr, "DEADSTEMC", fn, lno);

    fclose (bgc_file);
}

void ReadEPC (epctbl_struct *epctbl)
{
    int             i;
    char            fn[MAXSTRING];
    double          t1, t2, t3, t4, r1;
    FILE           *epc_file;
    char            cmdstr[MAXSTRING];

    epctbl->woody = (int *)malloc (NLCTYPE * sizeof (int));
    epctbl->evergreen = (int *)malloc (NLCTYPE * sizeof (int));
    epctbl->c3_flag = (int *)malloc (NLCTYPE * sizeof (int));
    epctbl->phenology_flag = (int *)malloc (NLCTYPE * sizeof (int));
    epctbl->onday = (int *)malloc (NLCTYPE * sizeof (int));
    epctbl->offday = (int *)malloc (NLCTYPE * sizeof (int));
    epctbl->transfer_days = (int *)malloc (NLCTYPE * sizeof (int));
    epctbl->litfall_days = (int *)malloc (NLCTYPE * sizeof (int));
    epctbl->leaf_turnover = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->froot_turnover = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->livewood_turnover = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->daily_mortality_turnover =
        (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->daily_fire_turnover =
        (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->alloc_frootc_leafc = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->alloc_newstemc_newleafc =
        (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->alloc_newlivewoodc_newwoodc =
        (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->alloc_crootc_stemc = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->alloc_prop_curgrowth =
        (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->avg_proj_sla = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->sla_ratio = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->lai_ratio = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->ext_coef = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->flnr = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->psi_open = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->psi_close = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->vpd_open = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->vpd_close = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->froot_cn = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->leaf_cn = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->livewood_cn = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->deadwood_cn = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->leaflitr_cn = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->leaflitr_flab = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->leaflitr_fucel = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->leaflitr_fscel = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->leaflitr_flig = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->frootlitr_flab = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->frootlitr_fucel = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->frootlitr_fscel = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->frootlitr_flig = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->deadwood_fucel = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->deadwood_fscel = (double *)malloc (NLCTYPE * sizeof (double));
    epctbl->deadwood_flig = (double *)malloc (NLCTYPE * sizeof (double));

    /* Read epc files */
    PIHMprintf (VL_VERBOSE, "\nRead ecophysiological constant files\n");

    for (i = 0; i < NLCTYPE; i++)
    {
        switch (i + 1)
        {
            case ENF:
                strcpy (fn, "input/epc/enf.epc");
                epc_file = fopen (fn, "r");
                break;
            case EBF:
                strcpy (fn, "input/epc/ebf.epc");
                epc_file = fopen (fn, "r");
                break;
            case DNF:
                strcpy (fn, "input/epc/dnf.epc");
                epc_file = fopen (fn, "r");
                break;
            case DBF:
                strcpy (fn, "input/epc/dbf.epc");
                epc_file = fopen (fn, "r");
                break;
            case GRASS:
                strcpy (fn, "input/epc/c3grass.epc");
                epc_file = fopen (fn, "r");
                break;
            case CLOSE_SHRUB:
                strcpy (fn, "input/epc/shrub.epc");
                epc_file = fopen (fn, "r");
                break;
            case OPEN_SHRUB:
                strcpy (fn, "input/epc/shrub.epc");
                epc_file = fopen (fn, "r");
                break;
            default:
                strcpy (fn, "N/A");
                epc_file = NULL;
        }

        if (strcasecmp (fn, "N/A") != 0)
        {
            CheckFile (epc_file, fn);
            PIHMprintf (VL_VERBOSE, " Reading %s\n", fn);

            /* Skip header file */
            fgets (cmdstr, MAXSTRING, epc_file);
            /* Read epc */
            /* woody/non-woody flag */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%d", &epctbl->woody[i]);
            /* evergreen/deciduous flag */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%d", &epctbl->evergreen[i]);
            /* C3/C4 flag */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%d", &epctbl->c3_flag[i]);
            /* transfer days */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%d", &epctbl->transfer_days[i]);
            /* litter fall days */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%d", &epctbl->litfall_days[i]);
            /* leaf turnover */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->leaf_turnover[i]);
            /* force leaf turnover fraction to 1.0 if deciduous */
            if (!epctbl->evergreen[i])
            {
                epctbl->leaf_turnover[i] = 1.0;
            }
            epctbl->froot_turnover[i] = epctbl->leaf_turnover[i];
            /* live wood turnover */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->livewood_turnover[i]);
            /* whole-plant mortality */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &t1);
            epctbl->daily_mortality_turnover[i] = t1 / 365.0;
            /* fire mortality */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &t1);
            epctbl->daily_fire_turnover[i] = t1 / 365.0;
            /* froot C:leaf C */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->alloc_frootc_leafc[i]);
            /* new stem C:new leaf C */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->alloc_newstemc_newleafc[i]);
            /* new livewood C:new wood C */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->alloc_newlivewoodc_newwoodc[i]);
            /* croot C:stem C */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->alloc_crootc_stemc[i]);
            /* new growth:storage growth */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->alloc_prop_curgrowth[i]);
            /* force storage growth to 0.0 if evergreen (following CLM-CN) */
            if (epctbl->evergreen[i])
            {
                epctbl->alloc_prop_curgrowth[i] = 1.0;
            }
            /* average leaf C:N */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->leaf_cn[i]);
            /* leaf litter C:N */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->leaflitr_cn[i]);
            /* test for leaflitter C:N > leaf C:N */
            if (epctbl->leaflitr_cn[i] < epctbl->leaf_cn[i])
            {
                PIHMprintf (VL_ERROR,
                    "Error: leaf litter C:N must be >= leaf C:N.\n");
                PIHMprintf (VL_ERROR,
                    "Change the values in epc file %s.\n", fn);
                PIHMexit (EXIT_FAILURE);
            }
            /* initial fine root C:N */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->froot_cn[i]);
            /* initial livewood C:N */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->livewood_cn[i]);
            /* initial deadwood C:N */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->deadwood_cn[i]);
            /* test for deadwood C:N > livewood C:N */
            if (epctbl->deadwood_cn[i] < epctbl->livewood_cn[i])
            {
                PIHMprintf (VL_ERROR,
                    "Error: livewood C:N must be >= deadwood C:N.\n");
                PIHMprintf (VL_ERROR,
                    "Change the values in epc file %s.\n", fn);
                PIHMexit (EXIT_FAILURE);
            }
            /* leaf litter labile proportion */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &t1);
            epctbl->leaflitr_flab[i] = t1;
            /* leaf litter cellulose proportion */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &t2);
            /* leaf litter lignin proportion */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &t3);
            epctbl->leaflitr_flig[i] = t3;

            /* test for litter fractions sum to 1.0 */
            if (fabs (t1 + t2 + t3 - 1.0) > FLT_COND_TOL)
            {
                PIHMprintf (VL_ERROR,
                    "Error: leaf litter proportions of labile, cellulose, "
                    "and lignin\n");
                PIHMprintf (VL_ERROR,
                    "must sum to 1.0. Check epc file %s and try again.\n",
                    fn);
                PIHMexit (EXIT_FAILURE);
            }
            /* calculate shielded and unshielded cellulose fraction */
            r1 = t3 / t2;
            if (r1 <= 0.45)
            {
                epctbl->leaflitr_fscel[i] = 0.0;
                epctbl->leaflitr_fucel[i] = t2;
            }
            else if (r1 > 0.45 && r1 < 0.7)
            {
                t4 = (r1 - 0.45) * 3.2;
                epctbl->leaflitr_fscel[i] = t4 * t2;
                epctbl->leaflitr_fucel[i] = (1.0 - t4) * t2;
            }
            else
            {
                epctbl->leaflitr_fscel[i] = 0.8 * t2;
                epctbl->leaflitr_fucel[i] = 0.2 * t2;
            }
            /* froot litter labile proportion */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &t1);
            epctbl->frootlitr_flab[i] = t1;
            /* froot litter cellulose proportion */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &t2);
            /* froot litter lignin proportion */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &t3);
            epctbl->frootlitr_flig[i] = t3;

            /* test for litter fractions sum to 1.0 */
            if (fabs (t1 + t2 + t3 - 1.0) > FLT_COND_TOL)
            {
                PIHMprintf (VL_ERROR,
                    "Error: froot litter proportions of labile, cellulose, "
                    "and lignin\n");
                PIHMprintf (VL_ERROR,
                    "must sum to 1.0. Check epc file %s and try again.\n",
                    fn);
                PIHMexit (EXIT_FAILURE);
            }
            /* calculate shielded and unshielded cellulose fraction */
            r1 = t3 / t2;
            if (r1 <= 0.45)
            {
                epctbl->frootlitr_fscel[i] = 0.0;
                epctbl->frootlitr_fucel[i] = t2;
            }
            else if (r1 > 0.45 && r1 < 0.7)
            {
                t4 = (r1 - 0.45) * 3.2;
                epctbl->frootlitr_fscel[i] = t4 * t2;
                epctbl->frootlitr_fucel[i] = (1.0 - t4) * t2;
            }
            else
            {
                epctbl->frootlitr_fscel[i] = 0.8 * t2;
                epctbl->frootlitr_fucel[i] = 0.2 * t2;
            }
            /* dead wood cellulose */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &t1);
            /* dead wood lignin */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &t2);
            epctbl->deadwood_flig[i] = t2;
            /* test for litter fractions sum to 1.0 */
            if (fabs (t1 + t2 - 1.0) > FLT_COND_TOL)
            {
                PIHMprintf (VL_ERROR,
                   "Error: deadwood proportions of cellulose and lignin "
                   "must sum\n");
                PIHMprintf (VL_ERROR,
                    "to 1.0. Check epc file %s and try again.\n", fn);
                PIHMexit (EXIT_FAILURE);
            }
            /* calculate shielded and unshielded cellulose fraction */
            r1 = t2 / t1;
            if (r1 <= 0.45)
            {
                epctbl->deadwood_fscel[i] = 0.0;
                epctbl->deadwood_fucel[i] = t1;
            }
            else if (r1 > 0.45 && r1 < 0.7)
            {
                t4 = (r1 - 0.45) * 3.2;
                epctbl->deadwood_fscel[i] = t4 * t1;
                epctbl->deadwood_fucel[i] = (1.0 - t4) * t1;
            }
            else
            {
                epctbl->deadwood_fscel[i] = 0.8 * t1;
                epctbl->deadwood_fucel[i] = 0.2 * t1;
            }
            /* canopy light ext coef */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->ext_coef[i]);
            /* all to projected LA ratio */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->lai_ratio[i]);
            /* canopy average projected specific leaf area */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->avg_proj_sla[i]);
            /* sunlit SLA ratio */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->sla_ratio[i]);
            /* Rubisco N fraction */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->flnr[i]);
            /* psi_sat */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->psi_open[i]);
            /* psi_close */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->psi_close[i]);
            /* vpd_max */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->vpd_open[i]);
            /* vpd_min */
            fgets (cmdstr, MAXSTRING, epc_file);
            sscanf (cmdstr, "%lf", &epctbl->vpd_close[i]);
        }
        else
        {
            epctbl->woody[i] = BADVAL;
            epctbl->evergreen[i] = BADVAL;
            epctbl->c3_flag[i] = BADVAL;
            epctbl->phenology_flag[i] = BADVAL;
            epctbl->onday[i] = BADVAL;
            epctbl->offday[i] = BADVAL;
            epctbl->transfer_days[i] = BADVAL;
            epctbl->litfall_days[i] = BADVAL;
            epctbl->leaf_turnover[i] = BADVAL;
            epctbl->froot_turnover[i] = BADVAL;
            epctbl->livewood_turnover[i] = BADVAL;
            epctbl->daily_mortality_turnover[i] = BADVAL;
            epctbl->daily_fire_turnover[i] = BADVAL;
            epctbl->alloc_frootc_leafc[i] = BADVAL;
            epctbl->alloc_newstemc_newleafc[i] = BADVAL;
            epctbl->alloc_newlivewoodc_newwoodc[i] = BADVAL;
            epctbl->alloc_crootc_stemc[i] = BADVAL;
            epctbl->alloc_prop_curgrowth[i] = BADVAL;
            epctbl->avg_proj_sla[i] = BADVAL;
            epctbl->sla_ratio[i] = BADVAL;
            epctbl->lai_ratio[i] = BADVAL;
            epctbl->ext_coef[i] = BADVAL;
            epctbl->flnr[i] = BADVAL;
            epctbl->psi_open[i] = BADVAL;
            epctbl->psi_close[i] = BADVAL;
            epctbl->vpd_open[i] = BADVAL;
            epctbl->vpd_close[i] = BADVAL;
            epctbl->froot_cn[i] = BADVAL;
            epctbl->leaf_cn[i] = BADVAL;
            epctbl->livewood_cn[i] = BADVAL;
            epctbl->deadwood_cn[i] = BADVAL;
            epctbl->leaflitr_cn[i] = BADVAL;
            epctbl->leaflitr_flab[i] = BADVAL;
            epctbl->leaflitr_fucel[i] = BADVAL;
            epctbl->leaflitr_fscel[i] = BADVAL;
            epctbl->leaflitr_flig[i] = BADVAL;
            epctbl->frootlitr_flab[i] = BADVAL;
            epctbl->frootlitr_fucel[i] = BADVAL;
            epctbl->frootlitr_fscel[i] = BADVAL;
            epctbl->frootlitr_flig[i] = BADVAL;
            epctbl->deadwood_fucel[i] = BADVAL;
            epctbl->deadwood_fscel[i] = BADVAL;
            epctbl->deadwood_flig[i] = BADVAL;
        }
    }
}

void ReadAnnFile (tsdata_struct *ts, char *fn)
{
    FILE           *fid;
    char            timestr[MAXSTRING];
    char            cmdstr[MAXSTRING];
    int             i;
    int             match;
    int             lno = 0;

    fid = fopen (fn, "r");
    CheckFile (fid, fn);
    PIHMprintf (VL_VERBOSE, " Reading %s\n", fn);

    ts->length = CountLine (fid, cmdstr, 1, "EOF");
    ts->ftime = (int *)malloc (ts->length * sizeof (int));
    ts->data = (double **)malloc (ts->length * sizeof (double *));

    FindLine (fid, "BOF", &lno, fn);
    for (i = 0; i < ts->length; i++)
    {
        ts->data[i] = (double *)malloc (sizeof (double));
        NextLine (fid, cmdstr, &lno);
        match =
            sscanf (cmdstr, "%s %lf", timestr, &ts->data[i][0]);

        if (match != 2)
        {
            PIHMprintf (VL_ERROR, "Error reading %s.\n", fn);
            PIHMprintf (VL_ERROR,
                "Please check file format near Line %s.\n", lno);
            PIHMexit (EXIT_FAILURE);
        }

        ts->ftime[i] = StrTime (timestr);
    }

    fclose (fid);
}

/*
 * restart_io.c
 * functions called to copy restart info between restart structure and
 * active structures
 * 
 * *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
 * Biome-BGC version 4.2 (final release)
 * See copyright.txt for Copyright information
 * *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
 */

#include "pihm.h"

void RestartInput (cstate_struct * cs, nstate_struct * ns, epvar_struct * epv, bgc_ic_struct * restart)
{
    cs->leafc = restart->leafc;
    cs->leafc_storage = restart->leafc_storage;
    cs->leafc_transfer = restart->leafc_transfer;
    cs->frootc = restart->frootc;
    cs->frootc_storage = restart->frootc_storage;
    cs->frootc_transfer = restart->frootc_transfer;
    cs->livestemc = restart->livestemc;
    cs->livestemc_storage = restart->livestemc_storage;
    cs->livestemc_transfer = restart->livestemc_transfer;
    cs->deadstemc = restart->deadstemc;
    cs->deadstemc_storage = restart->deadstemc_storage;
    cs->deadstemc_transfer = restart->deadstemc_transfer;
    cs->livecrootc = restart->livecrootc;
    cs->livecrootc_storage = restart->livecrootc_storage;
    cs->livecrootc_transfer = restart->livecrootc_transfer;
    cs->deadcrootc = restart->deadcrootc;
    cs->deadcrootc_storage = restart->deadcrootc_storage;
    cs->deadcrootc_transfer = restart->deadcrootc_transfer;
    cs->gresp_storage = restart->gresp_storage;
    cs->gresp_transfer = restart->gresp_transfer;
    cs->cwdc = restart->cwdc;
    cs->litr1c = restart->litr1c;
    cs->litr2c = restart->litr2c;
    cs->litr3c = restart->litr3c;
    cs->litr4c = restart->litr4c;
    cs->soil1c = restart->soil1c;
    cs->soil2c = restart->soil2c;
    cs->soil3c = restart->soil3c;
    cs->soil4c = restart->soil4c;
    cs->cpool = restart->cpool;

    ns->leafn = restart->leafn;
    ns->leafn_storage = restart->leafn_storage;
    ns->leafn_transfer = restart->leafn_transfer;
    ns->frootn = restart->frootn;
    ns->frootn_storage = restart->frootn_storage;
    ns->frootn_transfer = restart->frootn_transfer;
    ns->livestemn = restart->livestemn;
    ns->livestemn_storage = restart->livestemn_storage;
    ns->livestemn_transfer = restart->livestemn_transfer;
    ns->deadstemn = restart->deadstemn;
    ns->deadstemn_storage = restart->deadstemn_storage;
    ns->deadstemn_transfer = restart->deadstemn_transfer;
    ns->livecrootn = restart->livecrootn;
    ns->livecrootn_storage = restart->livecrootn_storage;
    ns->livecrootn_transfer = restart->livecrootn_transfer;
    ns->deadcrootn = restart->deadcrootn;
    ns->deadcrootn_storage = restart->deadcrootn_storage;
    ns->deadcrootn_transfer = restart->deadcrootn_transfer;
    ns->cwdn = restart->cwdn;
    ns->litr1n = restart->litr1n;
    ns->litr2n = restart->litr2n;
    ns->litr3n = restart->litr3n;
    ns->litr4n = restart->litr4n;
    ns->soil1n = restart->soil1n;
    ns->soil2n = restart->soil2n;
    ns->soil3n = restart->soil3n;
    ns->soil4n = restart->soil4n;
    ns->sminn = restart->sminn;
    ns->retransn = restart->retransn;
    ns->npool = restart->npool;

    epv->day_leafc_litfall_increment = restart->day_leafc_litfall_increment;
    epv->day_frootc_litfall_increment = restart->day_frootc_litfall_increment;
    epv->day_livestemc_turnover_increment = restart->day_livestemc_turnover_increment;
    epv->day_livecrootc_turnover_increment = restart->day_livecrootc_turnover_increment;
    epv->annmax_leafc = restart->annmax_leafc;
    epv->annmax_frootc = restart->annmax_frootc;
    epv->annmax_livestemc = restart->annmax_livestemc;
    epv->annmax_livecrootc = restart->annmax_livecrootc;
    epv->dsr = restart->dsr;

    epv->dormant_flag = restart->dormant_flag;
    epv->onset_flag = restart->onset_flag;
    epv->onset_counter = restart->onset_counter;
    epv->onset_gddflag = restart->onset_gddflag;
    epv->onset_fdd = restart->onset_fdd;
    epv->onset_gdd = restart->onset_gdd;
    epv->onset_swi = restart->onset_swi;
    epv->offset_flag = restart->offset_flag;
    epv->offset_counter = restart->offset_counter;
    epv->offset_fdd = restart->offset_fdd;
    epv->offset_swi = restart->offset_swi;
}

void RestartOutput (cstate_struct * cs, nstate_struct * ns, epvar_struct * epv, bgc_ic_struct *restart)
{
    restart->leafc = cs->leafc;
    restart->leafc_storage = cs->leafc_storage;
    restart->leafc_transfer = cs->leafc_transfer;
    restart->frootc = cs->frootc;
    restart->frootc_storage = cs->frootc_storage;
    restart->frootc_transfer = cs->frootc_transfer;
    restart->livestemc = cs->livestemc;
    restart->livestemc_storage = cs->livestemc_storage;
    restart->livestemc_transfer = cs->livestemc_transfer;
    restart->deadstemc = cs->deadstemc;
    restart->deadstemc_storage = cs->deadstemc_storage;
    restart->deadstemc_transfer = cs->deadstemc_transfer;
    restart->livecrootc = cs->livecrootc;
    restart->livecrootc_storage = cs->livecrootc_storage;
    restart->livecrootc_transfer = cs->livecrootc_transfer;
    restart->deadcrootc = cs->deadcrootc;
    restart->deadcrootc_storage = cs->deadcrootc_storage;
    restart->deadcrootc_transfer = cs->deadcrootc_transfer;
    restart->gresp_storage = cs->gresp_storage;
    restart->gresp_transfer = cs->gresp_transfer;
    restart->cwdc = cs->cwdc;
    restart->litr1c = cs->litr1c;
    restart->litr2c = cs->litr2c;
    restart->litr3c = cs->litr3c;
    restart->litr4c = cs->litr4c;
    restart->soil1c = cs->soil1c;
    restart->soil2c = cs->soil2c;
    restart->soil3c = cs->soil3c;
    restart->soil4c = cs->soil4c;
    restart->cpool = cs->cpool;

    restart->leafn = ns->leafn;
    restart->leafn_storage = ns->leafn_storage;
    restart->leafn_transfer = ns->leafn_transfer;
    restart->frootn = ns->frootn;
    restart->frootn_storage = ns->frootn_storage;
    restart->frootn_transfer = ns->frootn_transfer;
    restart->livestemn = ns->livestemn;
    restart->livestemn_storage = ns->livestemn_storage;
    restart->livestemn_transfer = ns->livestemn_transfer;
    restart->deadstemn = ns->deadstemn;
    restart->deadstemn_storage = ns->deadstemn_storage;
    restart->deadstemn_transfer = ns->deadstemn_transfer;
    restart->livecrootn = ns->livecrootn;
    restart->livecrootn_storage = ns->livecrootn_storage;
    restart->livecrootn_transfer = ns->livecrootn_transfer;
    restart->deadcrootn = ns->deadcrootn;
    restart->deadcrootn_storage = ns->deadcrootn_storage;
    restart->deadcrootn_transfer = ns->deadcrootn_transfer;
    restart->cwdn = ns->cwdn;
    restart->litr1n = ns->litr1n;
    restart->litr2n = ns->litr2n;
    restart->litr3n = ns->litr3n;
    restart->litr4n = ns->litr4n;
    restart->soil1n = ns->soil1n;
    restart->soil2n = ns->soil2n;
    restart->soil3n = ns->soil3n;
    restart->soil4n = ns->soil4n;
    restart->sminn = ns->sminn;
    restart->retransn = ns->retransn;
    restart->npool = ns->npool;

    restart->day_leafc_litfall_increment = epv->day_leafc_litfall_increment;
    restart->day_frootc_litfall_increment = epv->day_frootc_litfall_increment;
    restart->day_livestemc_turnover_increment = epv->day_livestemc_turnover_increment;
    restart->day_livecrootc_turnover_increment = epv->day_livecrootc_turnover_increment;
    restart->annmax_leafc = epv->annmax_leafc;
    restart->annmax_frootc = epv->annmax_frootc;
    restart->annmax_livestemc = epv->annmax_livestemc;
    restart->annmax_livecrootc = epv->annmax_livecrootc;
    restart->dsr = epv->dsr;

    restart->dormant_flag = epv->dormant_flag;
    restart->onset_flag = epv->onset_flag;
    restart->onset_counter = epv->onset_counter;
    restart->onset_gddflag = epv->onset_gddflag;
    restart->onset_fdd = epv->onset_fdd;
    restart->onset_gdd = epv->onset_gdd;
    restart->onset_swi = epv->onset_swi;
    restart->offset_flag = epv->offset_flag;
    restart->offset_counter = epv->offset_counter;
    restart->offset_fdd = epv->offset_fdd;
    restart->offset_swi = epv->offset_swi;
}

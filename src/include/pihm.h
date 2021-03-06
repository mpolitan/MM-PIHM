#ifndef PIHM_HEADER
#define PIHM_HEADER

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/stat.h>
#include <ctype.h>
#include <unistd.h>
#include <stdarg.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define VERSION     "0.5.0-alpha"

/*
 * SUNDIAL Header Files
 */
/* CVODE header file */
#include "cvode.h"
/* CVSPGMR linear header file */
#include "cvode_spgmr.h"
/* Definition of type N_Vector */
#ifdef _CVODE_OMP
#include "nvector_openmp.h"
#else
#include "nvector_serial.h"
#endif
/* UnitRoundoff, RSqrt, SQR functions */
#include "sundials_math.h"
/* CVDENSE header file */
#include "cvode_dense.h"

#ifdef _NOAH_
#include "spa.h"
#endif

#include "pihm_const.h"
#include "pihm_input_struct.h"
#include "elem_struct.h"
#include "river_struct.h"
#include "pihm_struct.h"
#include "pihm_func.h"
#endif

#include "pihm.h"

pihm_t_struct PIHMTime (int t)
{
    pihm_t_struct   pihm_time;
    struct tm      *timestamp;
    time_t          rawtime;

    rawtime = (time_t)t;
    timestamp = gmtime (&rawtime);

    pihm_time.t = t;
    pihm_time.year = timestamp->tm_year + 1900;
    pihm_time.month = timestamp->tm_mon + 1;
    pihm_time.day = timestamp->tm_mday;
    pihm_time.hour = timestamp->tm_hour;
    pihm_time.minute = timestamp->tm_min;
    strftime (pihm_time.str, 17, "%Y-%m-%d %H:%M", timestamp);

    return (pihm_time);
}

int StrTime (const char *timestr)
{
    struct tm      *timestamp;
    int             t;

    timestamp = (struct tm*) malloc (sizeof (struct tm));

    switch (strlen (timestr))
    {
        case 4:
            sscanf (timestr, "%d", &timestamp->tm_year);
            timestamp->tm_year -= 1900;
            timestamp->tm_mon = 0;
            timestamp->tm_mday = 1;
            timestamp->tm_hour = 0;
            timestamp->tm_min = 0;
            timestamp->tm_sec = 0;
            timestamp->tm_isdst = 0;
            timestamp->tm_gmtoff = 0;
            timestamp->tm_zone = 0;
            break;
        case 16:
            sscanf (timestr, "%d-%d-%d %d:%d", &timestamp->tm_year,
                &timestamp->tm_mon, &timestamp->tm_mday, &timestamp->tm_hour,
                &timestamp->tm_min);
            timestamp->tm_year -= 1900;
            timestamp->tm_mon--;
            timestamp->tm_sec = 0;
            timestamp->tm_isdst = 0;
            timestamp->tm_gmtoff = 0;
            timestamp->tm_zone = 0;
            break;
        default:
            PIHMprintf (VL_ERROR,
                "Error converting from time string to time.\n");
            PIHMexit (EXIT_FAILURE);
    }

    t = (int)timegm (timestamp);

    free (timestamp);

    return (t);
}




// collection of frequently used help functions
#include <time.h>
#include <ctype.h>

int time2tm(const int year, const int mon, const int mday, const int hour, const int min, const int sec, struct tm* tacq)
{
    setenv("TZ", "UTC", 1);
    tacq->tm_year = year - 1900;
    tacq->tm_mon  = mon - 1;
    tacq->tm_mday = mday;
    tacq->tm_hour = hour;
    tacq->tm_min  = min;
    tacq->tm_sec  = sec;
    time_t ttacq = mktime(tacq);
    tacq = gmtime(&ttacq);
    return tacq->tm_yday;
}

int date2tm(const int year, const int mon, const int mday, struct tm* tacq)
{
    return time2tm(year, mon, mday, 0, 0, 0, tacq);
}

double julian_day(const int year, const int month, const int mday, const int hour, const int min, const int sec, double msec)
{
    int jday = mday - 32075 + 1461 * (year + 4800 - (14 - month)/12)/4 +
	      367 * (month - 2 + (14 - month)/12*12)/12 -
	      3 * ((year + 4900 - (14 - month)/12)/100)/4;
    double Jday = (double) jday + (double) hour/24.0 + (double) min/1440.0 + ((double) sec + 0.001*msec)/86400;
    return Jday;
}

double jdn2jday(const int jday, const int hour, const int min, const int sec, double msec)
{
    double Jday = (double) jday + (double) hour/24.0 + (double) min/1440.0 + ((double) sec + 0.001*msec)/86400;
    return Jday;
}

short dayofyear(const int year, const int mon, const int mday)
{
    double yday = (julian_day(year, mon, mday, 0,0,0,0.0) - julian_day(year, 1, 1, 0,0,0,0.0) + 1.0);
    return  (short) yday;
}

int jday2tm(const double jday, double *msec, struct tm* tacq)
{
    int  julian_date, year, month, day, temp1, temp2;

    julian_date = (int) floor(jday);
    temp1 = julian_date + 68569;
    temp2 = 4*temp1/146097;
    temp1 = temp1 - (146097 * temp2 + 3) / 4;
    year  = 4000 * (temp1 + 1) / 1461001;
    temp1 = temp1 - 1461 * year/4 + 31;
    month = 80 * temp1 / 2447;
    day   = temp1 - 2447 * month / 80;
    temp1 = month / 11;
    month = month + 2 - 12 * temp1;
    year  = 100 * (temp2 - 49) + year + temp1;
    double dj = jday - (double) julian_date;
    int hour = (int) (dj * 24.0);
    int min  = ((int) (dj * 1440.0)) % 60;
    int sec  = ((int) (dj * 86400.0)) % 60;
    double dummy;
    if (msec) *msec = modf(dj * 86400.0, &dummy)*1000.0;
    return time2tm(year, month, day, hour, min, sec, tacq);
}

double dday(const short yday, const short hour, const short min, const short sec, const double msec)
{
    double day = (double) yday + hour/24.0 + min/1440.0 + (sec + 0.001*msec)/86400.0;
    return day;
}

int timestamp(const double yDay, const double offsetinsec, short* yday, short* hour, short* min, short* sec, short* msec, short* usec)
{
    double dummy;
    double day = yDay + offsetinsec/86400.0;
    *yday = (short) day;
    *hour = (short) (24.0*modf(day, &dummy));
    *min  = (short) (60.0*modf(24.0*day, &dummy));
    *sec  = (short) (60.0*modf(1440.0*day, &dummy));
    *msec = (short) (1000.0*modf(86400.0*day, &dummy));
    *usec = (short) (1000.0*modf(86400000.0*day, &dummy));
    return (int) (100000.0*day);
}

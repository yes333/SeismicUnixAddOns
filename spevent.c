/* Copyright (c) READ Well Services, 2010.*/
/* All rights reserved.                       */

/* SPEVENT: $Revision: 1.0 $ ; $Date: 2011/09/17 22:58:42 $	*/

#include "su.h"
#include "segy.h"
#include <segyhdr.h>
#include <cwp.h>
#include "header.h"
#include <signal.h>

/*********************** self documentation **********************/
char *sdoc[] = {
"                                                                               ",
" SPEVENT - multi-level micro-seismic event detection                           ",
"                                                                               ",
" spevent <stdin >stdout [optional parameters]                                  ",
"                                                                               ",
" Optional parameters:                                                          ",
"                                                                               ",
"   key=ep          gather sorting key                                          ",
"   nmax=12         max. number of receiver stations expected in a gather       ",
"   nc=3            number of components expected in a gather                   ",
"                                                                               ",
"   level=2         =1: compute background noise (rms)                          ",
"                   =2: detect, associate and separate (multiple) seismic events",
"                   =3: pick onset times of events (typically S-arrivals)       ",
"                   =4: determine S onset times of events based on SV/SH components",
"                                                                               ",
" Threshold regarding MEE amplitude for detecting event                         ",
"                                                                               ",
"   thres=3         rms amplitude ratio of three detection time window          ",
"   perc=75         percentage threshold for receivers with event detected      ",
"                                                                               ",
" Parameters to calculate modified energy envelop:                              ",
"                                                                               ",
"   win1=40  [ms]   window length above/before                                  ",
"   win2=20  [ms]   window length below/after                                   ",
"   win3=10  [ms]   window length below/after                                   ",
"   power1=1        power for rms ratio win2/win1                               ",
"   power3=1        power for rms ratio win3/trace                              ",
"                                                                               ",
" Parameters in addition for level 2 event correlation/output: 			",
"                                                                               ",
"   wina=10  [ms]   time difference event association/correlation window for adjacent receiver  ",
"   maxa=3*wina     max. event association/correlation window within same event ",
//"   ckey=counit     key holding the ordering no. of isolated events (1,2,3,...) ",
"                                                                               ",
"   wine=200 [ms]   trace length of output isolated seismic event               ",
"   rpos=0.6        relative position at output time window of correlated event ",
"                   0.5 means the center of output time window                  ",
"                                                                               ",
"   drop=0          =1: drop gather where no seismic event found                ",
//"                   =2: drop receiver where no clear seismic event found        ",
"   eeout=          compute and output modified energy envelop for event detection",
"                   =1 default for level 2                                      ",
"                   =0 default for level 3                                      ",
"                                                                               ",
" Parameters for level 3 S-arrival time picking:  				",
"                                                                               ",
"   t1=win1    [ms] begin of search window for S-arrivals                      ",
"   t2=ns-win2 [ms] end of search window for S-arrivals                        ",
"   tkey=laga       key holding S-arrival time in 1/10 ms                       ",
"                                                                               ",
"   verbose=0       >0  echo information                                        ",
"                                                                               ",
" Approach:                                                                     ",
" ---------                                                                     ",
" This program uses approach similar to first break detection (ShortTerm/LongTerm",
" ratio, i.e. RMS[win2/win1]), moreover, multiplied with a short term (win3) rms",
" over entire trace as weighting factor. To be more flexible, the two ratios are",
" powered by further parameters, creating a so called Modified Energy Envelop:  ",
"        MEE = RMS[win2/win2]^power1 * RMS[win3/trace]^power3                   ",
" The three windows, one before/above (win1) and two after/below (win2 & win3). ",
" If the MEE exceeds the threshold given by parameter (thres=), an event is said",
" to be found. All events found across all receivers over the whole panenl are  ",
" sorted by increasing arrival time and correlated to identify/isolate individual",
" seismic event (SE). Events of adjacent receivers are associated/correlated to ",
" belong to a same SE if their time differences are less than wina. The maximum ",
" time difference/gap allowed within same SE is specified by parameter (maxa=)  ",
" for any two consecutive events. Furthermore, a certain percentages (perc=) of ",
" receivers must have event detected to fully qualify a correlated event to be a",
" seismic event, which is then cut out and output with its median arrival time  ",
" positioned at the time specified by parameter rpos= within time window (wine=).",
" A unique event number ID is created by 100 * key + EventNo and stored to SU   ",
" keyword ep.                                                                   ",
" Input data must be sorted by receiver and component. If multiple SE are found ",
" during level=2 event detection the output must be sorted by keywords ep duse  ",
" before next level automatic arrival time picking.                             ",
" For level 1 seismic events are sought and rms of window length (wine=) around ",
" all events are subtracted from trace rms resulting in background noise rms for",
" individual receivers. These two rms values are saved to keywords fx and fy,   ",
" along with global background noise rms level over all receivers, computed and ",
" stored to keyword fz.                                                         ",
" Keyword corr holds the number of seismic events correlated.                   ",
"                                                                               ",
" Caveats:                                                                      ",
" the automatical picking of arrival times on level 3 requires relative neat MEE",
" local maximums within correlation window. It therefore fails more often if    ",
" multiple events or clear P and S arrivals overlap in time.                    ",
"                                                                               ",
" Example 1: background noise estimation                                        ",
"									        ",
" mirf2su length=0 file=myfile.rcd |\\",
" spevent key=fldr level=1 thres=2 |\\",
" suwind key=duse max=1 |\\",
" sugethw key=fldr,cdp,fx,fy,fz output=geom > rms-level.txt ",
"									        ",
" Example2: typical sequence for micro-seismic event detection/time pick        ",
"									        ",
" mirf2su length=8.1 overlap=0.1 file=$recno.rcd hfile=myfile.hdr |\\",
" spsort ep |\\",
" subfilt fstoplo=10 fpasslo=20 fpasshi=240 fstophi=360 |\\",
" sushw key=styp a=$StageNo |\\",
" sushw match=styp key=sx,sy,sdepth infile=injpoint.xyz |\\",
" sushw match=cdp key=gx,gy,gelev infile=receiver.xyz |\\",
" sushw match=cdp key=otrav infile=receiver-azimuth.tbl |\\",
" sp2crot key=cdp a=otrav mode=0 nc=3 pol=-1,1,-1 |\\",
" spevent key=ep level=2 nc=3 verbose=1 |\\",
" spsort ep |\\",
" spevent key=ep level=3 nc=4 verbose=1 |\\",
" tee $recno.su |\\",
" sugain pbal=1 |\\",
" spsort ep duse |\\",
" suxwigb perc=99.8 windowtitle=\"Record $recno after level 3 event time pick\" & ",
"                                                                               ",
" Version 1.2.0  last modified Nov. 2012 by Sanyu Ye				",
NULL
};

/* Credits:
 *
 * Sanyu Ye, sanyu_ye@yahoo.com
 *
 */
/**************** end self doc ***********************************/

// forward prototyp declaration
float calcEE(const int nc, const int nrcv, const int ns, float*** indata, float* rms, float** eedata);
void calcWinEE(const int nrcv, const int nw1, const int nw2, const int nw3, const int it1, const int it2,
    const float power1, const float power3, float* rms, float** eedata, float** eefunc);
float calcRms(const int itbeg, const int nt, const float* data);
int detectEvents(const int level, const int nrcv, const int ns,  const int nw1, const int nw2, const int nw3, const int nwa,
     const int it1, const int it2, const float threshold, float** evtfunc, int* nv, int** ite);

#include "sprinthlp.c"
#include "datetime.c"

int main(int argc, char **argv) 
{
    cwp_String key, ckey, tkey;          // header key word from segy.h
    cwp_String type, ctype, ttype;         // ... its type
    int index, cindex, tindex;          // ... its index
    Value val, valnew;                  // ... its value

    int level, totalEvents, hasEvent, nmax;

    int i, j, k;        // trace counter
    int nsegy;          // number of byte read
    int ntr, total, ngather, etotal; // number of traces and gathers etc
    int nsevent = 0;    // number of seismic events
    int nw1, nw2, nw3, nwa, nwe;       // number of samples
    int nrcv, ncmp;     // expected and actual number of receivers in gather
    int drop = 0;       // flag to output drop gather/receiver without event found
    int eeout = 0;      // flag to output trace of energy envelop
    int eof = 0;        // END OF File flag
    int verbose;        // flag for printing information

    float win1, win2, win3, wina, wine, maxa, rpos;   // time windows to calculate rms or average amplitudes		*/
    float dt;           // time sampling interval		*/
    float thres, power, power3; // amplitude ratio	*/
    float perc;         // percentage threshold of receivers with events detected
    float t1, t2, tevtm;
    float grms = 0.0;
    segy tr, outtr;

    /* Initialize */
    initargs(argc, argv);
    requestdoc(0);

    if (!getparint("verbose", &verbose)) verbose = 0;
    if (!getparint("nmax", &nmax))      nmax = 12;
    if (!getparint("nc", &ncmp))        ncmp = 3;
    if (!getparint("level", &level))   level = 2;
    if (!getparint("eeout", &eeout)) {
        eeout = (level > 2)? 0 : 1;
    }
    if (!getparint("drop", &drop))      drop = 0;
    if (!getparfloat("power1", &power))  power  = 1.0;
    if (!getparfloat("power3", &power3)) power3 = 1.0;
    if (!getparfloat("thres",  &thres)) thres =  3.0;
    if (!getparfloat("perc",  &perc))     perc =  75.0;
    if (!getparfloat("win1", &win1)) win1 = 40;
    if (!getparfloat("win2", &win2)) win2 = 20;
    if (!getparfloat("win3", &win3)) win3 = 10;
    if (!getparfloat("wina", &wina)) wina = 10;
    if (!getparfloat("maxa", &maxa)) maxa = 3*wina;
    if (!getparfloat("wine", &wine)) wine = 200;
    if (!getparfloat("rpos", &rpos)) rpos = 0.6;
    if (rpos < 0.2 || rpos > 0.8) {
        err("Center of correlated event should be place in the middle output window ( 0.2 < rpos < 0.8 )");
    }

    if(verbose > 1) {
        warn("level=%d  percentage=%2.0f  threshold=%3.1f   win1=%2.0f  win2=%2.0f  win3=%2.0f",
                level, perc, thres, win1, win2, win3);
        warn("correlation window=%2.0f  trace length=%3.0f", wina, wine);
    }

    /* Get info from first trace */
    if ((nsegy = gettr(&tr)) < HDRBYTES) err("can't read first trace");
    if (!tr.dt) err("dt header field must be set");
    dt = ((double) tr.dt) / 1000000.0;
    int ns = (int) tr.ns;

    /* Get sorting/gather key  */
    if (!getparstring("key", &key)) key = "ep";
    type = hdtype(key);
    index = getindex(key);
    gethval(&tr, index, &val);

    if (!getparstring("ckey", &ckey)) ckey = "counit";
    ctype = hdtype(ckey);
    cindex = getindex(ckey);
    if (!getparstring("tkey", &tkey)) tkey = "laga";
    ttype = hdtype(tkey);
    tindex = getindex(tkey);

    nw1 = NINT(0.001*win1/dt);
    nw2 = NINT(0.001*win2/dt);
    nw3 = NINT(0.001*win3/dt);
    nwa = NINT(0.001*wina/dt);
    nwe = NINT(0.001*wine/dt);
    const int nemax = ns / nw3 * 2;  // max. number of event
    if (!getparfloat("t1", &t1)) t1 = win1;
    if (!getparfloat("t2", &t2)) t2 = ns - win2 - 1;
    int it1 = (level > 2)? NINT(0.001 * t1 / dt) : nw1;
    int it2 = (level > 2)? NINT(0.001 * t2 / dt) : ns - MAX(nw2, nw3) - 1;
    if (it1 < nw1) it1 = nw1;
    if (it2 >= ns - MAX(nw2, nw3)) it2 = ns -  MAX(nw2, nw3) - 1;
    if (verbose>1 && level > 2) warn("search time window ms/ns  %1.0f/%d ~ %1.0f/%d", t1, it1, t2, it2);

    // allocate and reset memory for input and output traces
    segyhdr**  hdrs = (segyhdr**) ealloc2(ncmp, nmax, HDRBYTES);
    float*** indata = ealloc3float(ns, ncmp, nmax);
    float**  eedata = ealloc2float(ns, nmax);     // trace array for energy envelop
    float**   winee = ealloc2float(ns, nmax);     // trace array for energy envelop
    int**    ite = ealloc2int(nemax, nemax);          // event location in sample number
    int*      nv = ealloc1int(nmax);                   // number of event detected on every single receiver
    int*    hase = ealloc1int(nmax);           // flag for receiver with correlated event
    int*     nvn = ealloc1int(nemax);           // number of events of every seismic event
    int* itevent = ealloc1int(nemax);           // median time of events of a seismic event
    int*    ircv = ealloc1int(nemax*nemax);     // receiver where the event is detected
    int*    jevt = ealloc1int(nemax*nemax);     // serial event no within a receiver
    int*      ii = ealloc1int(nemax*nemax);     // serial index of events
    int**    iii = ealloc2int(nemax, nemax);    // index of event by seismic event no and receiver no
    float*     f = ealloc1float(nemax*nemax);   // time (ite) in float format of events
    float*  tevt = ealloc1float(nmax);          // cache for storing times of detected event at a receiver
    float*  rms = ealloc1float(nmax);          // cache for storing rms for every receiver over whole trace length
    float* bgrms = ealloc1float(nmax);         // cache for storing background rms for every receiver over whole trace length
    /* zero out data memory */
    memset(ii,   0, nemax * nemax * sizeof (int));
    memset(*iii, 0, nemax * nemax * sizeof (int));
    memset(ircv, 0, nemax * nemax * sizeof (int));
    memset(jevt, 0, nemax * nemax * sizeof (int));
    memset(*ite, 0, nemax * nemax * sizeof (int));
    memset(nv, 0, nmax * sizeof (int));
    memset(itevent, 0, nemax * sizeof (int));
    memset(*hdrs, 0, nmax * ncmp * HDRBYTES);
    memset(**indata, 0, nmax * ncmp * ns * FSIZE);
    memset( *eedata, 0, nmax * ns * FSIZE);
    memset(  *winee, 0, nmax * ns * FSIZE);
    memset(       f, 0, nemax * nemax * FSIZE);
    memset(     rms, 0, nmax * FSIZE);
    memset(   bgrms, 0, nmax * FSIZE);
    memset(    tevt, 0, nmax * FSIZE);
    
    /* Read headers and data while getting a count */
    ngather = ntr = total = etotal = 0;
    do {
        if (nsegy > HDRBYTES) gethval(&tr, index, &valnew);
        else eof = 1; //END_OF_FILE
        if (nsegy > HDRBYTES && !valcmp(type, val, valnew)) { /* same key and more data*/
            if (ntr > ncmp*nmax - 1) err("\nNumber of traces exceeding nmax=%d\n", nmax*ncmp);
            int ircv = ntr/ncmp;
            int icmp = ntr%ncmp;
            memcpy(&hdrs[ircv][icmp], &tr, HDRBYTES);
            memcpy(indata[ircv][icmp], tr.data, FSIZE * ns);
            ++ntr;
            val = valnew;
        } else { // new gather or END_OF_FILE
            int keyno = vtoi(type, val);

            nrcv = ntr/ncmp;

            ++ngather;
            if (verbose) warn("  processing %d traces %d-th gather (%s=%d)", ntr, ngather, key, keyno);

            if (ntr%ncmp !=0) {
                err(" Number of traces (%d != %d x %d) not equal number of receivers times compomnents", ntr, nrcv, ncmp);
            }
            
            total += ntr;

            if (level < 3 || ncmp <= 3 || level > 3) {
                grms = calcEE(ncmp, nrcv, ns, indata, rms, eedata); // calulate rms (3C) amplitude for every sample and receiver
                calcWinEE(nrcv, nw1, nw2, nw3, nw1, ns - nw2, power, power3, rms, eedata, winee);
            } else if (level == 3 && ncmp == 4) {
                 // copy from component 4
                if (eeout) {
                    for(i=0; i<nrcv; ++i) {
                        memcpy(eedata[i], indata[i][3], ns*FSIZE);
                        rms[i] = calcRms(0, ns, eedata[i]);
                    }
                    calcWinEE(nrcv, nw1, nw2, nw3, nw1, ns - nw2, power, power3, rms, eedata, winee);
                } else {
                    for(i=0; i<nrcv; ++i) memcpy(winee[i], indata[i][3], ns*FSIZE);
                }
            }
            
            totalEvents = detectEvents(level, nrcv, ns, nw1, nw2, nw3, nwa, it1 + 1, it2 - 1, thres, winee, nv, ite);
            if (verbose > 0) {
                fprintf(stderr, " Total %d events found for %d-th gather (%s=%d) with fldr=%d\n",
                        totalEvents, ngather, key, keyno, hdrs[0][0].fldr);
            }
            if (verbose > 1) {
                fprintf(stderr, "RCV NE  Time/RMS  Average RMS amplitude %5.3g in this gather\n", grms);
                for(i=0; i<nrcv; ++i) {
                    fprintf(stderr, "%2d  %2d", i+1, nv[i]);
                    for(j=0; j<nv[i]; ++j) fprintf(stderr, " %5.4f/%3.1f ", dt*ite[i][j], winee[i][ite[i][j]]);
                    fprintf(stderr, "\n");
                }
            }
            
            int nc;  // number of events correlated/belongs to previous event
            // sort event by event time
            for (i = 0, k = 0; i < nrcv; ++i) {
                for (j = 0; j < nv[i]; ++j) {
                    f[k] = (float) ite[i][j];   // event time in float
                    ircv[k] = i;                // rcv no event belongs to
                    jevt[k] = j;                // event no on that receiver
                    ii[k] = k;                  // ordering no
                    ++k;
                }
            }
            qkisort(k, f, ii);
            qksort(k, f);

            // figure out association/correlation
            nsevent = 0;
            int evtno = 0;
            for (i = 1; i < k;) {
                for (j = i, nc = 0, iii[nsevent][nc++] = evtno; j < k; ++j) {
                    i = j + 1;
                    ++evtno;
                    float maxdiff = MAX((float) nwa, MIN(maxa, ABS(ircv[ii[j]] - ircv[ii[j - 1]])*nwa));
                    if (f[j] - f[j - 1] <= maxdiff) {  // belong probably to same event
                        iii[nsevent][nc++] = evtno;
                    }
                    if (j == k - 1 || f[j] - f[j - 1] > maxdiff) { // last event or start of a new event
                        if (100.0 * nc < perc * nrcv) break;
                        int ie;
                        //for (ie=0; ie<nc; ++ie) { fprintf(stderr, " %d ", ircv[ii[iii[nsevent][ie]]] + 1); }
                        //fprintf(stderr, "\n");
                        memset(hase, 0, nmax * ISIZE);
                        for (ie = 0; ie < nc; ++ie) hase[ircv[ii[iii[nsevent][ie]]]] = 1;
                        int nrcvwe = 0;
                        for (ie = 0; ie < nrcv; ++ie) nrcvwe += hase[ie];
                        if (100.0 * nrcvwe >= perc * nrcv) { // an seismic event found
                            // print out info of seismic event
                            if (verbose > 1) {
                                fprintf(stderr, "\nSeismic Event No. %d by correlating %d small events (Time/EvtNo/RcvNo/EvtNoWithinReceiver)\n",
                                        nsevent + 1, nc);
                                for (ie = 0; ie < nc; ++ie) {
                                    int krcv = ircv[ii[iii[nsevent][ie]]];
                                    int kevt = jevt[ii[iii[nsevent][ie]]];
                                    //float tsort = f[iii[nsevent][ie]]*dt;
                                    float t = (float) ite[ircv[ii[iii[nsevent][ie]]]][jevt[ii[iii[nsevent][ie]]]]*dt;
                                    //if (tsort != t ) {
                                    //    warn("something wrong");
                                    //}
                                    fprintf(stderr, "%5.4f/%d/%d/%d ",
                                      t, ii[iii[nsevent][ie]], krcv + 1, kevt + 1);
                                }
                                fprintf(stderr, "\n");
                                //for (ie=0; ie<nc; ++ie) { fprintf(stderr, " %d/%d/%d ",iii[nsevent][ie], ii[iii[nsevent][ie]], ircv[ii[iii[nsevent][ie]]] + 1); }
                                //fprintf(stderr, "\n");
                            }

                            nvn[nsevent] = nc;
                            itevent[nsevent++] = (int) f[j - nc/2];  // median event time
                        }
                        break;
                    }
                }
            }
            etotal += nsevent;
            if(verbose > 0) {
                fprintf(stderr, "\n");
                fprintf(stderr, "   %d correlated seismic events found for gather (%s=%d)\n      with median time: ", nsevent, key, keyno);
                for(i = 0; i < nsevent; ++i) fprintf(stderr, "%5.3f ", itevent[i]*dt);
                fprintf(stderr, "\n");
            }

            // calculate background noises (trace rms - rms of event windows)
            float gbgrms = 0.0;
            if ( level == 1 ) {
                for (i=0; i<nrcv; ++i) {
                    if (nsevent == 0) {
                       bgrms[i] = rms[i];
                       continue;
                    }
                    int   nsig = 0;
                    float erms = 0;
                    for(j = 0; j < nsevent; ++j) {
                        int itb = MAX(0, itevent[j] - nwe/2);
                        int ite = MIN(ns - 1, itevent[j] + nwe/2);
                        int itn = ite - itb + 1;
                        float rmse = calcRms(itb, itn, eedata[i]);
                        erms += rmse*rmse*itn;
                        nsig += itn;
                    }
                    erms = rms[i]*rms[i]*ns - erms;
                    if (ns > nsig && erms > 0.0) bgrms[i] = sqrt(erms/(ns - nsig));
                    gbgrms += bgrms[i]*bgrms[i];
                }
                gbgrms = sqrtf(gbgrms/nrcv);
            }

            hasEvent = 0;
            int maxeno = 0;
            if (level > 2) {  // level = 3 or 4
                // try to pick that event which has most picks
                for (nc=0, i=0; i<nsevent; ++i) {
                    if (nvn[i] > nc ) {
                        maxeno = i;
                        nc = nvn[i];
                    }
                }
                tevtm = f[iii[maxeno][0] + nc*5/6]*dt;
                for (j=0, i=0; i<nc; ++i) {
                    //for (k=0; k<nc; ++k) { fprintf(stderr, " %d/%d/%d ",iii[maxeno][k], ii[iii[maxeno][k]], ircv[ii[iii[maxeno][k]]] + 1); }
                    //fprintf(stderr, "\n");
                    int krcv = ircv[ii[iii[maxeno][i]]];
                    float t = f[iii[maxeno][0] + i]*dt;
                    if (!tevt[krcv] || j < nrcv || t < tevtm) {
                        tevt[krcv] = t;
                        ++j;
                    } else {
                        int n = 0;
                        float tdiff1 = 0.0, tdiff2 = 0.0;
                        // check if it has smaller time tifferences to picks of its neiboring trace
                        if (krcv == 0) {  // first receiver
                            for (k=1; k<nrcv && n < 2; ++k) {
                                if (tevt[k]) {  // find non-null pick
                                    tdiff1 += (tevt[krcv] - tevt[k]);
                                    tdiff2 += (t - tevt[k]);
                                    ++n;
                                }
                            }
                            if ( tdiff1 < 0.0 || ABS(tdiff1) > ABS(tdiff2) ) {
                                tevt[krcv] = t;
                            }
                        } else if (krcv == nrcv-1) {  // last receiver
                            for (k=nrcv-2; k>=0 && n < 2; --k) {
                                if (tevt[k]) {  // find non-null pick
                                    tdiff1 += ABS(tevt[krcv] - tevt[k]);
                                    tdiff2 += ABS(t - tevt[k]);
                                    ++n;
                                }
                            }
                            if ( tdiff1 < 0.0 || ABS(tdiff1) > ABS(tdiff2) ) {
                                tevt[krcv] = t;
                            }
                       } else { // in middle
                            for (k=krcv-1; k>=0; --k) if (!tevt[k]) break;  // find non-null pick left
                            if (tevt[k]) {
                                tdiff1 += (tevt[krcv] - tevt[k]);
                                tdiff2 += (t          - tevt[k]);
                            }
                            for (k=krcv+1; k<nrcv; ++k) if (!tevt[k]) break;  // find non-null pick right
                            if (tevt[k]) {
                                tdiff1 += (tevt[krcv] - tevt[k]);
                                tdiff2 += (t          - tevt[k]);
                            }
                            if ( ABS(tdiff1) > ABS(tdiff2) ) {
                                tevt[krcv] = t;
                            }
                        }
                    }
                }
                if (verbose > 1) {
                    fprintf(stderr, "S-arrival time found:");
                    for(i=0; i<nrcv; ++i) fprintf(stderr, "%d/%5.4f  ", i+1, tevt[i]);
                    fprintf(stderr, "\n");
                }
                int n = 0;
                for(i=0; i<nrcv; ++i) if (tevt[i] > dt) ++n;
                if ( 100.0 * n >= nrcv * perc) hasEvent = 1;
            }

            // calc input trace starting time
            double tday = dday(hdrs[0][0].day, hdrs[0][0].hour, hdrs[0][0].minute, hdrs[0][0].sec, hdrs[0][0].timbas + 0.001*hdrs[0][0].trwf);
            if ( level >= 3 && drop == 1 && hasEvent == 0) {
                if(verbose) warn("Dropped %d-th gather (%s=%d) with fldr=%d",
                                 ngather, key, keyno, hdrs[0][0].fldr);
                continue;
            } else {
                int ncout = (eeout > 0)? ncmp + 1 : ncmp;
                for (i = 0; i < nrcv; ++i) {
                    for (j = 0; j < ncout; ++j) {
                        void* pdata = indata[i][j];
                        if (j == ncmp) { // ee component
                            pdata = winee[i];
                        }
                        if (j < ncmp) memcpy(&outtr, &hdrs[i][j], HDRBYTES);
                        memcpy(outtr.data, pdata, FSIZE * ns);
                        outtr.duse = j + 1;
                        outtr.corr = nsevent;
                        outtr.gaps = nv[i];
                        if (level == 1) {
                            outtr.unscale = gbgrms;
                            outtr.ungpow  = bgrms[i];
                            outtr.fz  = gbgrms;
                            outtr.fy  = bgrms[i];
                            outtr.fx  =   rms[i];
                            puttr(&outtr);
                        } else if (level == 2) {
                            for (k = 0; k < nsevent; ++k) {
                                int it = itevent[k] - NINT(nwe*rpos);  // put the event at related position of time window
                                int itstart = it;   // starting sample of input trace to be copied from
                                int ittr = 0;       // starting sample of output trace to be copied to
                                int nse = nwe;      // n of samples of output event trace
                                if (it < 0) {
                                    itstart = 0;
                                    ittr = -it;
                                    nse = nwe + it;
                                }
                                if ( it + nwe > ns ) nse -= it + nwe - ns; // if the end of window beyong trace end, cut short
                                double tte = 0.0;    // event time in 1/10 ms
                                int iv = 0;
                                for (iv = 0; iv < nvn[k]; ++iv) {
                                    if (i == ircv[ii[iii[k][iv]]]) {
                                        tte = 10000.0 * dt * (ite[i][jevt[ii[iii[k][iv]]]] - itstart);
                                    }
                                }
                                pdata = (j == ncmp)?  &winee[i][itstart] : &indata[i][j][itstart];
                                memset(outtr.data, 0, nwe*FSIZE);
                                memcpy(&outtr.data[ittr], pdata, nse * FSIZE);
                                outtr.ep = 100 * keyno + k + 1;  // uniq event number ID
                                outtr.ns = nwe;
                                // reset starting time
                                outtr.sdel = timestamp(tday, it*dt, &outtr.day, &outtr.hour, &outtr.minute,
                                              &outtr.sec, &outtr.timbas, &outtr.trwf);
                                setval(ctype, &val, (double) (k + 1));  // no. of isolated event
                                puthval(&outtr, cindex, &val);
                                setval(ttype, &val, tte);  // event time in 1/10 ms
                                puthval(&outtr, tindex, &val);
                                puttr(&outtr);
                            }
                        } else {  // level = 3 or 4
                            outtr.corr = hasEvent;
                            setval(ttype, &val, (double) 10000.0*tevt[i]);  // time of max event energy
                            puthval(&outtr, tindex, &val);
                            puttr(&outtr);
                        }
                    }
                }
            }
            val = valnew;
            // reset cache
            memset(ii,   0, nemax * nemax * sizeof (int));
            memset(*iii, 0, nemax * nemax * sizeof (int));
            memset(ircv, 0, nemax * nemax * sizeof (int));
            memset(jevt, 0, nemax * nemax * sizeof (int));
            memset(*ite, 0, nemax * nemax * sizeof (int));
            memset(nv, 0, nmax * sizeof (int));
            memset(itevent, 0, nemax * sizeof (int));
            memset(*hdrs, 0, nmax * ncmp * HDRBYTES);
            memset(**indata, 0, nmax * ncmp * ns * FSIZE);
            memset( *eedata, 0, nmax * ns * FSIZE);
            memset(  *winee, 0, nmax * ns * FSIZE);
            memset(       f, 0, nemax * nemax * FSIZE);
            memset(     rms, 0, nmax * FSIZE);
            memset(   bgrms, 0, nmax * FSIZE);
            memset(    tevt, 0, nmax * FSIZE);
            ntr = 0;
            continue;
        }
        nsegy = gettr(&tr);
    } while (!eof);

    if (verbose) {
        if (level == 1) warn("  %d events detected for total %d gathers of %d receivers", etotal, ngather, nrcv);
        else if (level == 2) warn("  %d correlated seismic events found for total %d gathers", etotal, ngather);
        else warn("Totally %d traces of %d gathers processed", total, ngather);
    }

    return (CWP_Exit());
}

float calcEE(const int nc, const int nrcv, const int ns, float*** indata, float* rms, float** eedata)
{
    int i, it;
    float sqr, rmspanel = 0.0;

    for(i=0; i<nrcv; ++i) {
        rms[i] = 0.0;
        for (it=0; it<ns; ++it) {
            if (nc < 5) {
                sqr = indata[i][0][it]*indata[i][0][it] +
                    indata[i][1][it]*indata[i][1][it] + indata[i][2][it]*indata[i][2][it];
            } else {
                sqr = indata[i][nc-2][it]*indata[i][nc-2][it] + indata[i][nc-1][it]*indata[i][nc-1][it];
            }
            eedata[i][it] = sqrtf(sqr);
            rms[i] += sqr;
        }
        rmspanel += rms[i];
        rms[i] = sqrtf(rms[i]/ns);
    }
    return sqrtf(rmspanel/(nrcv*ns));
}

void calcWinEE(const int nrcv, const int nw1, const int nw2, const int nw3, const int it1, const int it2,
        const float power, const float power3, float* rms,  float** eedata, float** evtfunc)
{
    int i, it;
    float sqr1, sqr2, sqr3;

    for(i=0; i<nrcv; ++i) {
        sqr3 = sqr2 = sqr1 = 0.0;
        for (it=it1; it<=it2; ++it) {
            sqr1 = calcRms(it - nw1, nw1, eedata[i]);
            sqr2 = calcRms(it, nw2, eedata[i]);
            sqr3 = calcRms(it, nw3, eedata[i]);
            evtfunc[i][it] = (sqr1 == 0.0)? 0.0 : powf(sqr3/rms[i], power3)*powf(sqr2/sqr1, power);
        }
    }
    return;
}

int detectEvents(const int level, const int nrcv, const int ns,  const int nw1, const int nw2, const int nw3, const int nwa,
     const int it1, const int it2, const float threshold, float** evtfunc, int* nv, int** ite)
{
    int i, j, it, ntotal = 0;

    for(i=0; i<nrcv; ++i) {
        j = 0;
        for (it=it1; it<=it2; ++it) {
            if (evtfunc[i][it] > threshold) {
                if (evtfunc[i][it - 1] > evtfunc[i][it] || evtfunc[i][it + 1] > evtfunc[i][it]) {  // must on local maximum
                    continue;
                }
                if(j > 0 && it - ite[i][j-1] < nw2) {
                    if ( it - ite[i][j-1] < nwa && evtfunc[i][it] > evtfunc[i][ite[i][j-1]] ) {
                        ite[i][j-1] = it;  // replace previous pick
                        it += 2;
                        continue;
                    } else if ( evtfunc[i][it] < 0.5*evtfunc[i][ite[i][j-1]] ) {  // ignore
                        it += 2;
                        continue;
                    }
                }
                ++nv[i];
                ite[i][j++] = it;
                it += 2;  // skip to next window
            }
        }
        if (nv[i] == 0 && level > 2) {
            // try to find the maximum instead
            float eemax = 0.0;
            for (it=it1; it<=it2; ++it) {
                if ( evtfunc[i][it] > eemax ) {
                    eemax = evtfunc[i][it];
                    ite[i][0] = it;
                }
            }
            nv[i] = 1;
        }
        ntotal += nv[i];
    }

    return ntotal;
}

float calcRms(const int itbeg, const int nt, const float* data)
{
    float sum = 0.0;
    int nc, it;
    for (it = 0, nc = 0; it < nt; ++it ) {
        if (data[itbeg + it] != 0.0 ) {
            ++nc;
            sum += (data[itbeg + it])*(data[itbeg + it]);
        }
    }
    return (nc > 0)? sqrtf(sum/nc) : 0.0;
}


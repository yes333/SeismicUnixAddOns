

sedis2sac verbose=1 fimg=sedis.img rcvno=27 length=8 



sedis2sac verbose=1 fimg=sedis.img rcvno=27 length=8 


	
proj +proj=utm  +lon_0=113e +ellps=WGS84 -r -E <<EOF
11.50092 113.53701  27
 9.17118 114.36068  40
EOF

	
proj +proj=utm  +lon_0=113e +ellps=WGS84 -E <<EOF
113.53701 11.50092  27
 114.36068  9.17118 40
EOF


ln -sf /local/data/OBS973-1/obs27/prof1.obs27p1.img input.img

sedis2sac rcvno=27 fimg=input.img lat=11.50092 long=113.53701 \
         gy=1272589.59 gx=776758.06 t0=2009,4,15 length=72 verbose=1


ln -sf /local/data/OBS973-1/obs36/prof1.obs36p1.img input.img

sedis2sac rcvno=36 fimg=input.img \
         t0=2009,4,15 length=96 verbose=1


ln -sf /local/data/OBS973-1/obs40/prof1.obs40p1.img input.img

sedis2sac rcvno=40 fimg=input.img lat=9.17118 long=114.36068 \
         gy=1015505.70 gx=869412.91 t0=2009,4,15 length=72 verbose=1


ln -s /local/data/OBS973-1/ukooa1.new1 ukooa.txt 

awk '{yday=substr($0, 5, 3); shot=substr($0, 20, 5); \
      lat=substr($0, 26, 2) + substr($0, 28, 2)/60. + substr($0, 30, 5)/3600.; \
      lon=substr($0, 36, 3) + substr($0, 39, 2)/60. + substr($0, 41, 5)/3600.; \
      printf("%7.5f %8.5f 0 %4d 2009 04 %02d %02d %02d %6.3f\n", \
      lat, lon, shot, yday-90, substr($0, 8, 2), substr($0, 10, 2), substr($0, 12, 6));}' \
      ukooa.txt |\
proj +proj=utm  +lon_0=113e +ellps=WGS84 -r |\
awk '{printf("S %4d %.1f %.1f %.1f %4d %02d %02d %02d %02d %6.3f\n", \
     $4, $1, $2, $3, $5, $6, $7, $8, $9, $10)}' |\
cat > shot.tbl

#more shot.tbl

head shot.tbl

tail shot.tbl

#

ln -sf ../sedis/shot.tbl nav.tbl

for stn in 27 36 40
do

for cmp in 1 2 3 4 
do

out=obs$stn-$cmp.su


ln -sf `ls ../sedis/obs$stn-*.sac$cmp` data.sac

# t0=0.116 correction for digital filter delay

sac2su fsac=data.sac fnav=nav.tbl t0=0.116 verbose=4 rv=6 |\
sushw key=f1 a=0 > $out


suwind key=offset max=140000  < $out |\
sufilter f=2,4,20,30 |\
sushw key=unscale match=offset values=0,5,5000,5,200000,200 interp=1 |\
suhtmath op=mult key=unscale |\
#sugain pbal=1 |\
suchw key1=gdel key2=sy key3=gy c=-1 |\
sushw key=ungpow match=gdel interp=1 values=-1,-1,1,1 |\
suchw key1=offset key3=offset key4=ungpow b=0 c=1 |\
suchw key1=f2 key2=offset d=1000 |\
suxwigb key=f2 perc=90 wbox=1400 hbox=880 grid1=solid \
    windowtitle="OBS $stn  Component $cmp" &


suwind key=offset max=200000  < $out |\
sufilter f=2,6,30,60 |\
sushw key=unscale match=offset values=200000,200,5000,5,0,5 interp=1 |\
suhtmath op=mult key=unscale |\
sugain pbal=1 |\
suximage perc=98 wbox=1400 hbox=880 grid1=solid \
    windowtitle="OBS $stn  Component $cmp" &

done

done



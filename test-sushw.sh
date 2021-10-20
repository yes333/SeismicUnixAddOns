
suwind key=cdpt min=0 max=0 accept=1001,1040,1071,1100,1201,1300,1400,1431,1462,1477,1478 < 2h.su | suwind key=nhs min=0 max=0 accept=1,2,5,20,35,50,69,75,79,80 > ! test.su

suwind key=cdpt min=0 max=0 accept=1001,1040,1071,1100,1201,1300,1400,1431,1462,1477,1478 < ../sptest/2hzxy-ch120.su | suwind key=nhs min=0 max=0 accept=1,2,5,20,35,50,69,75,79,80 | susort2 ep nhs cdpt > 3c.su




/home/yes/apps/read/sushw/sushw  match=cdpt,nhs key=gdel,sdel values=1400,5,1400,5,1400,75,1400,75,1100,5,1100,5,1100,75,1100,75 \
	sort=-1,1 interp=1 verbose=111 < test.su |\
sugethw key=cdpt,nhs,gdel,sdel |\
more


/home/yes/apps/read/sushw/sushw  match=cdpt,nhs key=gdel,sdel values=1100,5,1100,5,1100,75,1100,75,1400,5,1400,5,1400,75,1400,75 \
	sort=1,1 interp=1 verbose=111 < test.su |\
sugethw key=cdpt,nhs,gdel,sdel |\
more

/home/yes/apps/read/sushw/sushw  match=cdpt,nhs key=gdel,sdel values=1100,75,1100,75,1100,5,1100,5,1400,75,1400,75,1400,5,1400,5 \
	sort=1,-1 interp=1  < test.su |\
sugethw key=cdpt,nhs,gdel,sdel |\
more

/home/yes/apps/read/sushw/sushw  match=nhs,cdpt key=gdel,sdel values=75,1150,1150,75,75,1400,1400,75,5,1100,1100,5,5,1400,1400,5 \
	sort=-1,1 interp=1  < test.su |\
sugethw key=cdpt,nhs,gdel,sdel |\
more


/home/yes/apps/read/sushw/sushw  match=nhs,cdpt key=gdel,sdel values=75,1400,1400,75,75,1100,1100,75,5,1400,1400,5,5,1100,1100,5 \
	sort=-1,-1 interp=1 verbose=1 < test.su |\
sugethw key=cdpt,nhs,gdel,sdel |\
more


/home/yes/apps/read/sushw/sushw  match=ep,nhs,cdpt key=gdel,sdel \
 values=1,75,1400,1400,75,1,75,1100,1100,75,1,5,1400,1400,5,1,5,1100,1100,5,5,75,1400,7000,375,5,75,1100,5500,375,5,5,1400,7000,25,5,5,1100,5500,25 \
	sort=1,-1,-1 interp=1 verbose=10 < 3c.su |\
sugethw key=ep,nhs,cdpt,gdel,sdel |\
more


function [imgsv,imgsh] = HFM_migrate_3d(data,tstep,freqlow,freqhgh,rcsr,xi,yi,zi,search_area,vpfact,tp,p_ray,vsvfact,tsv,sv_ray,vshfact,tsh,sh_ray)
%
% [imgsv,imgsh] = HFM_migrate_3d(data,tstep,freqlow,freqhgh,rcsr,xi,yi,zi,search_area,vpfact,tp,p_ray,vsvfact,tsv,sv_ray,vshfact,tsh,sh_ray);
%
% migrate/deconvolve 3C data by reconstructing the compressional and shear
% source signal at each point in image space and perform a
% semblance-weighted deconvolution of one with respect to the other.
%
% The 3d migration is done in a "rolodex" fashion, in radial sections:
% finding the direction to each of the receivers and project the data to  
% and perpendicular to the plane from the source point through each of the 
% receivers. The elastic wave field is assumed to have traveled in this 
% plane, with P and Sv in the plane, and  Sh polarized perpendicular to the 
% plane. This will happen if the velocities are rotationally symmetric. 
% In this case the travel times are independent of the azimuth, and the 
% index into the pre-calculated 2D (xz) travel-time matrix is determined by
% the horizontal distance travelled. If the receiver array is not vertical, 
% this index would be different for each receiver. 
%
% To save on processing time, the event is assumed to be within a distance 
% ("search area") from the injection point. The xyz coordinates of the 
% injection point is contained in the first 3 elements of the array rcsr 
% (which has 6 X number of receivers)  
%
%
% INPUT variables
% data        - data (frequency domain)
%               dimensions(samples,receiver_levels,3)
%
% tstep       - sampling interval (in sec)
% freqlow     - low frequency  (in Hz)
% freqhgh     - high frequency (in Hz)
%
% rcsr        - xyz source (used of injection point coordinates), xyz receivers, 
%               dimensions(6,receiver_levels) 
%
% xi          - row vector of x (East)     values for image
% yi          - row vector of y (North)    values for image
% zi          - row vector of z (Vertical) values for image
%
% search_area - only generate the image within a sphere of radius "search
%               area" around the presumed injection point
%
% vpfact      - scaling factor divided into the travel times for compressional signal
% tp          - precomputed travel times for compressional signal in given velocity model
%               dimensions(length(zi),length(xi),size(data,2))
% p_ray       - precomputed ray angles at receiver for compressional signal
%               dimensions(length(zi),length(xi))
%
% vsvfact     - scaling factor divided into the travel times for vertical shear signal
% tsv         - precomputed travel times for vertical shear signal in given velocity model
%               dimensions(length(zi),length(xi))
% sv_ray      - precomputed ray angles at receiver for vertical shear signal
%               dimensions(length(zi),length(xi))
%
% vshfact     - scaling factor divided into the travel times for horizontal shear signal
% tsh         - precomputed travel times for horizontal shear signal in given velocity model
%               dimensions(length(zi),length(xi))
% sh_ray      - precomputed ray angles at receiver for horizontal shear signal
%               dimensions(length(zi),length(xi))
%
% OUTPUT variables
% imgsv       - vertical shear image
%               dimensions(length(zi),length(xi),length(yi))
% imgsh       - horizontal shear image
%               dimensions(length(zi),length(xi),length(yi))
%
% Jakob Haldorsen, August 27, 2012
%

[ltr,ntr,nchan] = size(data);

domega = -2*pi/ltr;
lom  = [0:ltr-1];
lom(floor(ltr/2+1):ltr) = lom(floor(ltr/2+1):ltr) - ltr ;
omega  = domega*lom';

% low and high frequency relative to Nyquist
flow  = 2*freqlow*tstep;
fhigh = 2*freqhgh*tstep;

% low and high sample to use
fsample = (floor(flow*ltr/2)+1);
lsample = min([floor(ltr/2) (floor(fhigh*ltr/2)+1)]);

% pull out the coordinates for the injection point
source = rcsr(1:3,1);

% allocate space for a subset of time and ray parameters
nzi = length(zi);
tpi    = zeros(nzi,ntr,1);
tvi    = zeros(nzi,ntr,1);
thi    = zeros(nzi,ntr,1);
pi_ray = zeros(nzi,ntr,1);
vi_ray = zeros(nzi,ntr,1);
hi_ray = zeros(nzi,ntr,1);

% traveltimes, velocity-scaled, fractional samples
tp  =  tp/(tstep*vpfact);
tsv = tsv/(tstep*vsvfact);
tsh = tsh/(tstep*vshfact);

% distances from injection point
distx2 = (xi-source(1)).^2;
disty2 = (yi-source(2)).^2;
distz2 = (zi-source(3)).^2;
distx  = sqrt(distx2);
disty  = sqrt(disty2);
distz  = sqrt(distz2);

% find the volume that could possibly contain a frac' event 
% distance from image point to receiver location
% these are the range of indices for the appropriate subset of image space
iiix = find(distx<=search_area);
iiiy = find(disty<=search_area);
iiiz = find(distz<=search_area);

% allocate image space for separat Sv and Sh images 
imgsv = zeros(length(zi),length(xi),length(yi));
imgsh = zeros(length(zi),length(xi),length(yi));

% Fourrier transform data
%scr   = fft(data);
scr   = data;

% select subset of data for imaging
scr   = scr(fsample:lsample,:,:);
omega = omega(fsample:lsample);
scr2  = scr;
ltr = size(scr2,1);
for ix=iiix
    for iy=iiiy;
        distxy = distx2(ix)+disty2(iy);
        if sqrt(distxy)<=search_area
            xt  = [xi(ix)-rcsr(4,:,1);yi(iy)-rcsr(5,:,1)];
            xt2 = sqrt(sum(xt.^2,1));
            if xt2~=0
                for ii=1:ntr
                    cs = xt(1,ii)./xt2(ii);
                    sn = xt(2,ii)./xt2(ii);
                    rotc = [cs sn;-sn cs];
                    scr2(:,ii,[2 3]) = squeeze(scr(:,ii,[2 3]))*rotc';
                    xt3 = rotc*xt(:,ii);
                    [a,iix] = min(abs(xt3(1)-xi));
                    rix     = max(1,min(size(p_ray,2),iix));
                    tpi(:,ii)  =       tp(:,rix,ii);
                    tvi(:,ii)  =      tsv(:,rix,ii);
                    thi(:,ii)  =      tsh(:,rix,ii);
                    pi_ray(:,ii) =  p_ray(:,rix,ii);
                    vi_ray(:,ii) = sv_ray(:,rix,ii);
                    hi_ray(:,ii) = sh_ray(:,rix,ii);
                end
            else
                for ii=1:ntr
                    [a,iix] = min(abs(xi-0));
                    rix     = max(1,min(size(p_ray,2),iix));
                    tpi(:,ii)  =       tp(:,rix,ii);
                    tvi(:,ii)  =      tsv(:,rix,ii);
                    thi(:,ii)  =      tsh(:,rix,ii);
                    pi_ray(:,ii) =  p_ray(:,rix,ii);
                    vi_ray(:,ii) = sv_ray(:,rix,ii);
                    hi_ray(:,ii) = sh_ray(:,rix,ii);
                end
            end
            % "rix" is the index into the travel-time matrices.
  
            sr  = rcsr(6,:,1);
            sc1 = scr2(:,:,1);
            sc2 = scr2(:,:,2);
            sc3 = scr2(:,:,3);
            scrv = zeros(nzi,1);
            scrh = zeros(nzi,1);

            for iz=iiiz;%1:nzi
                dists = sqrt(distz2(iz)+distxy);
                if dists<=search_area;
                    dist  = sqrt((zi(iz)-sr).^2+xt2.^2);
   
                    t0      = min(tpi(iz,:));
                    
                    phasep  = exp(-1i*omega*(tpi(iz,:)-t0));
                    phasesv = exp(-1i*omega*(tvi(iz,:)-t0));
                    phasesh = exp(-1i*omega*(thi(iz,:)-t0));

                    scrp  = (sc2.*(ones(ltr,1)*(dist.*(cos(pi_ray(iz,:))))) + sc1.*(ones(ltr,1)*(dist.*(sin(pi_ray(iz,:)))))).*phasep;            
                    scrsv = (sc2.*(ones(ltr,1)*(dist.*(sin(vi_ray(iz,:))))) - sc1.*(ones(ltr,1)*(dist.*(cos(vi_ray(iz,:)))))).*phasesv;
                    scrsh = (sc3.*(ones(ltr,1)*(dist))).*phasesh;           

                    eng   = mean(conj(scrp ).*scrp ,2);
                    engv  = mean(conj(scrsv).*scrsv,2);
                    engh  = mean(conj(scrsh).*scrsh,2);

                    scrp  = median(real(scrp ),2)-1i*median(imag(scrp ),2);%mean(scrp,2);%
                    scrsv = median(real(scrsv),2)+1i*median(imag(scrsv),2);%mean(scrsv,2);%
                    scrsh = median(real(scrsh),2)+1i*median(imag(scrsh),2);%mean(scrsh,2);%

                    scrsv = scrsv.*scrp./(eng+engv+engh);
                    scrsv = (2*real(mean(scrsv))).^2;
                    %scrsv = (2*real(mean(scrsv(fsample:lsample)))).^2;
                    scrv(iz) = scrsv(1);

                    scrsh = scrsh.*scrp./(eng+engv+engh);
                    scrsh = (2*real(mean(scrsh))).^2;
                    %scrsh = (2*real(mean(scrsh(fsample:lsample)))).^2;
                    scrh(iz) = scrsh(1);
                else
                    scrv(iz) = 0;
                    scrh(iz) = 0;
                end
            end
            imgsv(:,ix,iy) = scrv;
            imgsh(:,ix,iy) = scrh;
        end
    end
end



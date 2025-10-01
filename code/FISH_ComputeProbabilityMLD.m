

function [mld, ctrl] = FISH_ComputeProbabilityMLD(mld, fish, sst)

% Gather information on seasonal productivity and MLD (inshore/offshore)

% Seasonal cycle of CCS productive phase
prod =[-1.2 -1.2 -1.1 -.5 0 0.6 1.2 1.4 1.3 0.7 -.1 -.85];
prod = prod-mean(prod); prod= prod/std(prod);
ctrl.prod = prod;

% Estiamted ranges from close inspection of data 
mm = [2 4.5 7 9.5 11.5 14 16.5 19 21.5 23.5];
mldmin = [45 14 10 12 30 45 14 10 12 30];
mldmax = [120 44 20 27 65 120 44 20 27 65];
mldmin2 = interp1(mm, mldmin, [2:13]);
mldmax2 = interp1(mm, mldmax, [2:13]);
mldmin2 = [mldmin2(end) , mldmin2(1:end-1)];
mldmax2 = [mldmax2(end) , mldmax2(1:end-1)];
ctrl.mld_mins = mldmin2;
ctrl.mld_maxs = mldmax2;

mu=meanNaN(fish.ta.data,1); sigma=stdNaN(fish.ta.data,1)*1.2;
mint=mu-sigma*2;
maxt=mu+sigma*2;

for imon =1:12    
    % mask the region outside fish temperatur range
    sstmean = sst.seas(:,:,imon);
    sstmean2 = interp2(sst.lon', sst.lat', sstmean', mld.lon, mld.lat);
    sstmean2(sstmean2 > maxt) = nan;
    sstmean2(sstmean2 < mint) = nan;
    sstmean2(~isnan(sstmean2))=1;    
    % mask the region outside fish mld range
    tmp = mld.seas(:,:,imon);
    tmp(tmp>fish.mld.maxs(imon)) =nan;
    tmp(tmp<fish.mld.mins(imon)) =nan;
    % combine the two masks
    tmp=tmp.*sstmean2;
    % estimate mld depth for inshore and offshore within the masked area
    in = find(mld.lon > -200 & mld.lon < -150);
    mld.offshore(imon)= meanNaN(tmp(in),1);
    in = find(mld.lon > -140);
    mld.inshore(imon)= meanNaN(tmp(in),1);     
end

% add mix layer cycle to seasonal control array
ctrl.mld_offshore = mld.offshore;
ctrl.mld_inshore = mld.inshore;

% find esges of CCS productive system
sstmean = interp2(sst.lon', sst.lat', meanNaN(sst.seas(:,:,5:9),3)', mld.lon, mld.lat);
mldmean = meanNaN(mld.seas(:,:,5:9),3);

in = find(mldmean < 20 & sstmean <maxt & sstmean > mint);
tmp = mldmean*nan;
tmp(in) = mldmean(in);
clf;
FISH_2dmap(tmp, mld);
set(gca,'ylim',[15 62])
for i=1:12
    p = fish2d(i);
    month = str2num(
    in = find (p.m



PLOTTING =0;
if PLOTTING == 1
% MLD gradient, direction, and amplitude
[Gx, Gy, grad_angle, amplitude] = computeGradient(mld.mean);
[I,J] = size(mld.lon);
dx = 40; I=1:dx:I; J=1:dx:J;

figure(1)
clf
subplot(2,1,1)
plot_2dmap(mld.mean, mld)
title ('Mean MLD 1980-2024')
set(gca, 'FontSize',14)
subplot(2,1,2)
plot_2dmap(grad_angle, mld)
quiver(mld.lon(I,J),mld.lat(I,J),Gx(I,J),Gy(I,J),4,'k')
colormap("parula")
title 'Direction of MLD Gradient (angle 90 degree = North)'
set(gca, 'FontSize',14)
end

% compute gradient as a function of month
for imon =1:12
    [Gx(:,:,imon), Gy(:,:,imon), mld.grad_angle(:,:,imon), ~] = ...
        computeGradient(mld.seas(:,:,imon));
%     clf
%     subplot(2,1,1)
%     plot_2dmap(mld.seas(:,:,imon), mld)
%     title (imon)
%     subplot(2,1,2)
%     plot_2dmap(mld.grad_angle(:,:,imon), mld)
%     quiver(mld.lon(I,J),mld.lat(I,J),Gx(I,J,imon),Gy(I,J,imon),4,'k')
%     gradsmap4
%     pause
end

TEMP = mld.seas;
PROB = mld.seas;
PROB(:)=1; mu=meanNaN(fish.ta.data,1); sigma=stdNaN(fish.ta.data,1)*1.2;
in = find (TEMP(:) < mu+sigma*1 & TEMP(:) > mu-stdNaN(fish.ta.data,1)*1); 
PROB(in) = 100;
in = find (TEMP(:) >= mu+sigma*1); PROB(in) = -60;
in = find (TEMP(:) >= mu+sigma*2); PROB(in) = -30;
in = find (TEMP(:) >= mu+sigma*3); PROB(in) = -10;
in = find (TEMP(:) >= mu+sigma*4); PROB(in) = -1;
in = find (TEMP <= mu-sigma*1); PROB(in) = 60;
in = find (TEMP <= mu-sigma*2); PROB(in) = 30;
in = find (TEMP <= mu-sigma*3); PROB(in) = 10;
in = find (TEMP(:) <= mu-sigma*4); PROB(in) = 1;

mld.prob = PROB;
in = find(isnan(mld.seas));
mld.prob(in)=1;
in = find(isnan(mld.grad_angle));
mld.grad_angle(in)=180;

prob = mld.mask;
prob(~isnan(prob))=100;
prob(isnan(prob))=1;
mld.prob_land = prob;





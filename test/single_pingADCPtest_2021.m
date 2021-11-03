%Use toolbox to process single ping ADCP data.
% problem with current version of toolbox is that the plotting of
% single-ping data has a huge time overhead. Goal is to use the toolbox as
% much as possible to get the data processed, but avoid the plotting until
% the end - or all together.
%
% Bec Cowley, May/June, 2020
%Have updated workhorseParse to do beam to ENU conversion
% 2021: adapted to review screening of single-ping data with fish detection
% and correlation magnitude

% Set up:
clear
  % get the toolbox execution mode
  mode = readProperty('toolbox.mode');
toolboxVersion = ['2.6.13 - ' computer];
auto = false;
iMooring = [];
nonUTCRawData = importManager(toolboxVersion, auto, iMooring);
%% first do mapping and an ENU conversion so we can see velocities in the screening plots
%run using the pp window to ensure flags are handled properly;
% run these routines:
% adcpBinMappingPP
% adcpWorkhorseCorrMagPP - RDI recommended threshold = 64
% adcpWorkhorseEchoRangePP - RDI recommended threshold = 50
% adcpWorkhorseVelocityBeam2EnuPP
for k = 1:length(nonUTCRawData), nonUTCRawData{k}.meta.index = k; end
% preprocess data
 [sden, cancel] = preprocessManager(nonUTCRawData, 'qc', mode, false); 
% magneticDeclinationPP
%% Determine thresholds using some low-overhead plots:

%first check the single-ping screening thresholds (echo range and
%corrMagPP)
% If these are not suitable, need to adjust in the pp step (next cell).
% Might need to delete the *.ppp files in the raw_data folder before
% running PP as the toolbox will use whatever is in the ppp file, regardless
% of what you type into the starting box.

adcpScreeningThresholds(sden,1)
return
%% Now re-do with ensemble averaging; apply pp tools in this order:
% adcpBinMappingPP
% adcpWorkhorseCorrMagPP - RDI recommended threshold = 64
% adcpWorkhorseEchoRangePP - RDI recommended threshold = 50
% adcpWorkhorseVelocityBeam2EnuPP
% magneticDeclinationPP
% adcpWorkhorseRDIensembleAveragingPP - adjust the averaging interval,
%           default = 60 minutes
% depthPP

% preprocess data
 [ppData, cancel] = preprocessManager(nonUTCRawData, 'qc', mode, false); 
%  [nnData, cancel] = preprocessManager(nonUTCRawData, 'qc', mode, false); 

%%
% now look at ensemble-averaged data for appropriate autoQC thresholds to
% be applied in the autoQCManager step.
%error velocity histograms etc
%cmag already screened. Only apply if not done in single ping screening
%above, combined with surface test

% adcpEnsemblesThresholds(ppData,inst_index,ervthreshold,hvel,vvel,tilt,echo amp,cmagthreshold)
adcpEnsemblesThresholds(ppData,3,0.06,1,0.05,40,[],[])

return
%% Run autoQC tests, apply selected thresholds to see impact:
% do one at a time so that the results plot shows the fails for each
% individual test.
% run QC routines over raw data
aqc = autoQCManager(ppData);

% otherwise return new QC data
autoQCData = aqc;

% now review the results of the flags applied with each test
% comment out as required
adcpThresholdsResults(autoQCData,'surfacetest')
% adcpThresholdsResults(autoQCData,'vvel')
% adcpThresholdsResults(autoQCData,'echorange')
% adcpThresholdsResults(autoQCData,'erv')
% adcpThresholdsResults(autoQCData,'cmag')
% adcpThresholdsResults(autoQCData,'echo')
return
%% now we have the thresholds, run the toolbox as is, with the values determined.
% need to delete the pqc file for the RDIs before doing the imos toobox as
% I think that the values in pqc are overwriting edits in the
% configuration. need to check. YES, BUG/

%Next step, export to single ping netcdf and then re-import to matlab, bin
%average and stack.

%% checking what is going on with the bin mapping
% use echo intensity
nn = 1;
tr = nonUTCRawData{nn}.dimensions{1}.data;
ipresrel = getVar(nonUTCRawData{nn}.variables, 'PRES_REL');
pressure = nonUTCRawData{nn}.variables{ipresrel}.data;
iheight = getVar(nonUTCRawData{nn}.dimensions, 'DIST_ALONG_BEAMS');
Bins    = nonUTCRawData{nn}.dimensions{iheight}.data';
depr = repmat(-1*gsw_z_from_p(pressure,-27),1,length(Bins)) - repmat(Bins,length(tr),1);
for a = 1:4
    eval(['ech = nonUTCRawData{nn}.variables{' num2str(a+8) '}.data;'])
    eval(['echp = ppData{nn}.variables{' num2str(a+8) '}.data;'])
    eval(['echn = nnData{nn}.variables{' num2str(a+8) '}.data;'])
    df = diff(ech,[],2);
    dfp = diff(echp,[],2);
    dfn = diff(echn,[],2);
    figure(2);clf
    subplot(311)
    pcolor(tr,depr(:,2:end)',df');shading flat;axis ij;grid;caxis([0 40])
    title(['Echo Raw ' num2str(a)])
    subplot(312)
    pcolor(tr,depr(:,2:end)',dfp');shading flat;axis ij;grid;caxis([0 40])
    title(['Echo Toolbox map ' num2str(a)])
    subplot(313)
    pcolor(tr,depr(:,2:end)',dfn');shading flat;axis ij;grid;caxis([0 40])
    title(['Echo New map ' num2str(a)])
    linkaxes
    print(['Echocomp' num2str(a) '_down.png'],'-dpng')
%     pcolor(tr,depr',echp'-echn');shading flat;axis ij;grid;caxis([-0.5 0.5])
    pause
end

%% comparison of ENU conversions with velocity conversion to enu
% need to have processed 12474 (down) and 16531 (up) with toolbox convert
% to ENU, no bin mapping. In nnData structure. {1} is 12474 and {2} is
% 16531
% variables are 'names', 'info', 'wt','sens'
% nn = 1;
% load('/oa-decadal-climate/work/observations/oceanobs_data/EACdata/mooring/EAC1909_2105/raw_data/RDI/12474000.000.mat')
nn = 2;
load('/oa-decadal-climate/work/observations/oceanobs_data/EACdata/mooring/EAC1909_2105/raw_data/RDI/16431000.000.mat')
% match the times
ti = nnData{nn}.dimensions{1}.data;
tiv = sens.time/24/60/60 + datenum('1970/01/01','yyyy/mm/dd/');

ii = find(ti >= tiv(1)-5/60/60/24 & ti <= tiv(end)+5/60/60/24);

%now compare the v and u and w
u = nnData{nn}.variables{5}.data(ii,:);
v = nnData{nn}.variables{6}.data(ii,:);
w = nnData{nn}.variables{7}.data(ii,:);
uv = squeeze(wt.vel(:,:,1));
vv = squeeze(wt.vel(:,:,2));
wv = squeeze(wt.vel(:,:,3));

%display the differences
range(u-uv)
range(v-vv)
range(w-wv)

%plot
dist = nnData{nn}.dimensions{2}.data;
ipresrel = getVar(nnData{nn}.variables, 'PRES_REL');
pressure= nnData{nn}.variables{ipresrel}.data(ii);
depth = repmat(-1*gsw_z_from_p(pressure,-27),1,length(dist)) - repmat(dist,1,length(ii))';
figure(1);clf
pcolor(tiv,depth',(v-vv)');axis ij;shading flat;grid;colorbar

% Bugger all difference for both instruments, but not identical. Had to change
% orientation of 12474 in velocity by adding 180deg to roll. Velocity
% didn't recognise that it is a down-facing instrument. But both have
% differences. Let's accept it and move on.
%% comparison of ENU conversions +cell mapping with velocity conversion to enu
% need to have processed 12474 (down) and 16531 (up) with toolbox convert
% to ENU, no bin mapping. In nnData/ppData structures, {1} is 12474 and {2} is
% 16531
% ppData has the original toolbox bin mapping method, nnData is the new
% toolbox bin mapping method
% variables are 'names', 'info', 'wt','sens'
nn = 1;
load('/oa-decadal-climate/work/observations/oceanobs_data/EACdata/mooring/EAC1909_2105/raw_data/RDI/12474000.000map.mat')
% nn = 2;
% load('/oa-decadal-climate/work/observations/oceanobs_data/EACdata/mooring/EAC1909_2105/raw_data/RDI/16431000.000map.mat')
% match the times
ti = ppData{nn}.dimensions{1}.data;
tiv = sens.time/24/60/60 + datenum('1970/01/01','yyyy/mm/dd/');

ii = find(ti >= tiv(1)-5/60/60/24 & ti <= tiv(end)+5/60/60/24);

%now compare the v and u and w
u = ppData{nn}.variables{5}.data(ii,:);
v = ppData{nn}.variables{6}.data(ii,:);
w = ppData{nn}.variables{7}.data(ii,:);
un = nnData{nn}.variables{5}.data(ii,:);
vn = nnData{nn}.variables{6}.data(ii,:);
wn = nnData{nn}.variables{7}.data(ii,:);
uv = squeeze(wt.vel(:,:,1));
vv = squeeze(wt.vel(:,:,2));
wv = squeeze(wt.vel(:,:,3));

%display the differences
range(u-uv)
range(v-vv)
range(w-wv)

%plot
dist = ppData{nn}.dimensions{2}.data;
ipresrel = getVar(ppData{nn}.variables, 'PRES_REL');
pressure= ppData{nn}.variables{ipresrel}.data(ii);
depth = repmat(-1*gsw_z_from_p(pressure,-27),1,length(dist)) - repmat(dist,1,length(ii))';
figure(1);clf
subplot(211)
pcolor(tiv,depth',(u-un)');axis ij;shading flat;grid;colorbar;caxis([-0.5 0.5])
title('Original - Velocity')
subplot(212)
pcolor(tiv,depth',(un-uv)');axis ij;shading flat;grid;colorbar;caxis([-0.5 0.5])
title('New - Velocity')
linkaxes

figure(2);clf
subplot(211)
pcolor(tiv,depth',uv');axis ij;shading flat;grid;colorbar;caxis([-0.5 0.5])
title('Velocity')
subplot(212)
pcolor(tiv,depth',un');axis ij;shading flat;grid;colorbar;caxis([-0.5 0.5])
title('New')
linkaxes

figure(3);clf
subplot(121)
vdif = u-uv;
histogram(vdif(:),[-3:0.02:4])
hold on
vdif = un-uv;
histogram(vdif(:),[-3:0.02:4])
grid
title('U Original - Velocity')
subplot(122)
vdif = v-vv;
histogram(vdif(:),[-3:0.02:4])
hold on
vdif = vn-vv;
histogram(vdif(:),[-3:0.02:4])
grid
title('V Original - Velocity')
linkaxes

% very different in places
if nn == 1 %compare downward and aquadopp
    ti1 = nnData{1}.dimensions{1}.data;
    u = ppData{nn}.variables{5}.data;
    v = ppData{nn}.variables{6}.data;
    un = nnData{nn}.variables{5}.data;
    vn = nnData{nn}.variables{6}.data;
    figure(4);clf
    subplot(211);hold on
    plot(ti,v(:,23))
    plot(ti,vn(:,23))
    plot(ti1,nnData{1}.variables{5}.data,'k')
    title('V')
    subplot(212);hold on
    plot(ti,u(:,23))
    plot(ti,un(:,23))
    plot(ti1,nnData{1}.variables{6}.data,'k')
    title('U')
    linkaxes
    legend('Old','New','AqD')
end

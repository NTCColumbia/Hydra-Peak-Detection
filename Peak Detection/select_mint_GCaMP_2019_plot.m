% ANALYZE MEAN INTENSITY TRACES: RP/CB WINDOW SELECTION, THRESHOLD
% DEFINITION AND EVENT DETECTION

format compact
%% set parameters
clear;
fpath = 'C:\Users\rylab\Desktop\Codes for Osmo paper\Peak Detection\'; % path to the intensity file 
fname = '06_Osm_iGCaMP_H03_ContHM_06261700006'; % name of the intensity txt file %exclude text
bkgname = ''; % '' if no bkg file
spath = fpath; % path to save result
figpath = fpath;
framerate = 2; %fps

%% process data
% read data
mint = dlmread([fpath fname '.txt'],'\t',1,1);
mint = reshape(mint,1,[]);
if ~isempty(bkgname)
    bkg = reshape(dlmread([fpath bkgname '.txt'],'\t',1,1),1,[]);
    mint = mint-bkg;
end


% dF/F
% dff = calc_dff(mint);
dff = calc_dff_mw(mint,200);

% filter
% dff_filt = fir_lowpass_ct(dff,0.3); % 0.15
% dff_filt = fir_lowpass_ct(dff,0.15); % 0.15
dff_filt = dff;

% first derivative
ddff = gradient(dff_filt);
numT = length(dff_filt);

figure; set(gcf,'color','w')
subplot(2,1,1); plot(dff_filt); ylabel('filtered dF/F'); xlim([1 numT])
subplot(2,1,2); plot(ddff); ylabel('d(dF/F)'); xlabel('Time(Min)'); xlim([1 numT])

%% define time bins
% if 0
ysc = [min(dff_filt) max(dff_filt)];
% ysc = [0.8,1.2];

% select cb time bin
figure;  set(gcf,'color','w','position',[200 200 1534 500])
plot(dff_filt)
xlim([1 numT]);ylim([ysc(1) ysc(2)])
xlabel('Time(Min)');ylabel('dF/F')
title('select CB time bin; right click to stop')
% [xx,yy] = ginput;
% cb_time_bin = reshape(xx,[],2);
cb_time_bin = select_timebin(gcf);
cb_time_bin(cb_time_bin(:,1)<=0,1) = 1;
cb_time_bin(cb_time_bin(:,2)>numT,2) = numT;
close(gcf)

% select rp time bin
figure; set(gcf,'color','w','position',[200 200 1534 500])
hold on
for ii = 1:size(cb_time_bin,1)
    patch([cb_time_bin(ii,1),cb_time_bin(ii,2),cb_time_bin(ii,2),cb_time_bin(ii,1),cb_time_bin(ii,1)],...
        [ysc(1) ysc(1) ysc(2) ysc(2) ysc(1)],0.8*[1 1 1],'edgecolor','none','facealpha',0.7);
end
plot(dff_filt);xlim([1 numT]); ylim([ysc(1) ysc(2)])
xlabel('Frame');ylabel('d(dF/F)')
title('select RP time bin; right click to stop')
% [xx,yy] = ginput;
% rp_time_bin = reshape(xx,[],2);
rp_time_bin = select_timebin(gcf);
close(gcf)

% if selected time point is not within recording time
cb_time_bin(cb_time_bin>length(dff)) = length(dff);
cb_time_bin(cb_time_bin<=0) = 1;
rp_time_bin(rp_time_bin>length(dff)) = length(dff);
rp_time_bin(rp_time_bin<=0) = 1;

% display time bins
figure; set(gcf,'color','w','position',[200 200 1534 500])
hold on
for ii = 1:size(cb_time_bin,1)
    h1 = patch([cb_time_bin(ii,1),cb_time_bin(ii,2),cb_time_bin(ii,2),cb_time_bin(ii,1),cb_time_bin(ii,1)],...
        [ysc(1) ysc(1) ysc(2) ysc(2) ysc(1)],0.8*[1 1 1],'edgecolor','none','facealpha',0.7);
end
for ii = 1:size(rp_time_bin,1)
    h2 = patch([rp_time_bin(ii,1),rp_time_bin(ii,2),rp_time_bin(ii,2),rp_time_bin(ii,1),rp_time_bin(ii,1)],...
        [ysc(1) ysc(1) ysc(2) ysc(2) ysc(1)],[1 0.8 0.8],'edgecolor','none','facealpha',0.7);
end
plot(dff_filt,'k'); xlim([1 numT]); ylim([ysc(1) ysc(2)])
xlabel('Frame');ylabel('dF/F')
legend([h1 h2],'CB','RP')
saveas(gcf,[figpath fname '_time_bin.fig']);

% end
%keyboard;

%% CB peaks
ysc = [min(ddff) max(ddff)];
% ysc=1*[-0.01,0.04];

% select threshold
figure; set(gcf,'color','w','position',[200 200 1534 500])
hold on
for ii = 1:size(cb_time_bin,1)
    patch([cb_time_bin(ii,1),cb_time_bin(ii,2),cb_time_bin(ii,2),cb_time_bin(ii,1),cb_time_bin(ii,1)],...
        [ysc(1) ysc(1) ysc(2) ysc(2) ysc(1)],0.8*[1 1 1],'edgecolor','none','facealpha',0.7);
end
plot(ddff);xlim([1 numT]);ylim([ysc(1) ysc(2)])
xlabel('Frame');ylabel('d(dF/F)')
title('select CB threshold')
thresh_cb = select_thresh(gcf);
thresh_cb = sort(thresh_cb,'ascend');
close(gcf)

% detect peaks
cb_locs = [];
for ii = 1:size(cb_time_bin,1)
    [pks,clocs] = findpeaks(ddff(cb_time_bin(ii,1):cb_time_bin(ii,2)));
    keepindx = pks>thresh_cb(1) & pks<thresh_cb(2);
    clocs = clocs(keepindx)+cb_time_bin(ii,1);
    cb_locs(end+1:end+length(clocs)) = clocs;
end

figure; set(gcf,'color','w','position',[200 200 1534 500])
hold on
for ii = 1:size(cb_time_bin,1)
    patch([cb_time_bin(ii,1),cb_time_bin(ii,2),cb_time_bin(ii,2),cb_time_bin(ii,1),cb_time_bin(ii,1)],...
        [ysc(1) ysc(1) ysc(2) ysc(2) ysc(1)],0.8*[1 1 1],'edgecolor','none','facealpha',0.7);
end
plot(ddff);
plot([1 numT],thresh_cb(1)*[1 1],'k--');
plot([1 numT],thresh_cb(2)*[1 1],'k--');
xlabel('Frame');ylabel('d(dF/F)')
xlim([1 numT]);ylim([ysc(1) ysc(2)])

%% RP peaks
figure; set(gcf,'color','w','position',[200 200 1534 500])
hold on
for ii = 1:size(rp_time_bin,1)
    patch([rp_time_bin(ii,1),rp_time_bin(ii,2),rp_time_bin(ii,2),rp_time_bin(ii,1),rp_time_bin(ii,1)],...
        [ysc(1) ysc(1) ysc(2) ysc(2) ysc(1)],0.8*[1 1 1],'edgecolor','none','facealpha',0.7);
end
plot(ddff);xlim([1 numT]);ylim([ysc(1) ysc(2)])
xlabel('Frame');ylabel('d(dF/F)')
title('select RP threshold')
thresh_rp = select_thresh(gcf);
thresh_rp = sort(thresh_rp,'ascend');
close(gcf)

% detect peaks
rp_locs = [];
for ii = 1:size(rp_time_bin,1)
    [pks,clocs] = findpeaks(ddff(rp_time_bin(ii,1):rp_time_bin(ii,2)));
    keepindx = pks>thresh_rp(1) & pks<thresh_rp(2);
    clocs = clocs(keepindx)+rp_time_bin(ii,1);
    rp_locs(end+1:end+length(clocs)) = clocs;
end

% make time vector and look up locs of events
endtime = (1/framerate) * (length(mint)-1) ./ 60;
t = linspace(0,endtime,length(mint));
cb_locs_time = t(cb_locs);
rp_locs_time = t(rp_locs);
cb_time_bin = t(cb_time_bin);
rp_time_bin = t(rp_time_bin);

%% plot result
figure; set(gcf,'color','w','position',[255 336 1510 647])
% cb
ysc_filt = [min(dff_filt) max(dff_filt)];
subplot(2,2,1); hold on
for ii = 1:size(cb_time_bin,1)
    patch([cb_time_bin(ii,1),cb_time_bin(ii,2),cb_time_bin(ii,2),cb_time_bin(ii,1),cb_time_bin(ii,1)],...
        [ysc_filt(1) ysc_filt(1) ysc_filt(2) ysc_filt(2) ysc_filt(1)],0.8*[1 1 1],'edgecolor','none','facealpha',0.7);
end
plot(t, dff_filt);
scatter(cb_locs_time,dff_filt(cb_locs),'r.');
xlim([1 endtime]);
ylim(ysc_filt)
ylabel('dF/F'); title('CB')
subplot(2,2,3); hold on
% ysc = [min(ddff) max(ddff)];
for ii = 1:size(cb_time_bin,1)
    patch([cb_time_bin(ii,1),cb_time_bin(ii,2),cb_time_bin(ii,2),cb_time_bin(ii,1),cb_time_bin(ii,1)],...
        [ysc(1) ysc(1) ysc(2) ysc(2) ysc(1)],0.8*[1 1 1],'edgecolor','none','facealpha',0.7);
end
plot(t, ddff);
scatter(cb_locs_time,ddff(cb_locs-1),'r.');
xlim([1 endtime]);
ylim(ysc)
xlabel('Time (min)')
ylabel('d(dF/F)')
% rp
subplot(2,2,2); hold on
for ii = 1:size(rp_time_bin,1)
    patch([rp_time_bin(ii,1),rp_time_bin(ii,2),rp_time_bin(ii,2),rp_time_bin(ii,1),rp_time_bin(ii,1)],...
        [ysc_filt(1) ysc_filt(1) ysc_filt(2) ysc_filt(2) ysc_filt(1)],[1 0.8 0.8],'edgecolor','none','facealpha',0.7);
end
plot(t, dff_filt);
scatter(rp_locs_time,dff_filt(rp_locs),'r.');
xlim([1 endtime]);
ylim(ysc_filt)
ylabel('dF/F'); title('RP')
subplot(2,2,4); hold on
for ii = 1:size(rp_time_bin,1)
    patch([rp_time_bin(ii,1),rp_time_bin(ii,2),rp_time_bin(ii,2),rp_time_bin(ii,1),rp_time_bin(ii,1)],...
        [ysc(1) ysc(1) ysc(2) ysc(2) ysc(1)],[1 0.8 0.8],'edgecolor','none','facealpha',0.7);
end
plot(t, ddff);
scatter(rp_locs_time,ddff(rp_locs-1),'r.');
xlim([1 endtime]);
ylim(ysc)
xlabel('Time (min)')                   %JS: convert time to min
ylabel('d(dF/F)')

saveas(gcf,[figpath fname '_detection.fig']);

%% display result
% ldisp(table((1:length(rp_locs))',rp_locs','variablenames',{'number','peak_time'}))

if 0
%% add rp2
rp2_locs = []; % [vector], or [], modify manually
rp2_locs = intersect(rp_locs,rp2_locs);
rp1_locs = setdiff(rp_locs,rp2_locs);
rp2_manual = []; % put the undetected rp2 here
rp2_locs = [rp2_locs,rp2_manual]; % final result is a concatenation of the automated ones and manual ones

%% Make Raster plots

Raster_for_CBRP_073017

saveas(gcf,[figpath fname '_raster.fig']);

%% Calculate CB RP total peaks, duration, frequency 

cb_peaks = length(cb_locs_time)   % #of CB peaks

cb_length1 = size(cb_time_bin); cb_length = cb_length1(1);

cb_total_length=0;
for i = 1:cb_length
   cb_calc = cb_time_bin(i,2) - cb_time_bin(i,1);
   cb_total_length = cb_total_length + cb_calc; 
end
cb_total_length          % total duration of CB session

cb_frequency = cb_peaks / cb_total_length   % CB firing frequency per minute

rp_peaks = length(rp_locs_time)   % #of rp peaks

rp_length1 = size(rp_time_bin); rp_length = rp_length1(1);

rp_total_length=0;
for i = 1:rp_length
   rp_calc = rp_time_bin(i,2) - rp_time_bin(i,1);
   rp_total_length = rp_total_length + rp_calc; 
end
rp_total_length          % total duration of rp session

rp_frequency = rp_peaks / rp_total_length   % rp firing frequency per minute

%% sliding window for RP raster plot practice



duration = 120; % in minutes for the recording 
windowsize = 2; %minutes - 15 to match chris     % play around here
framerate = 5; %fps                                          THIS PARAMETER CAN BE CHANGED
frame_number =  120*60*framerate;     %number of total frame, can input this %initially movie param number

endtime = (1/framerate) * (frame_number-1) ./ 60;  %create time vector
t = linspace(0,endtime,frame_number);

windowframes = framerate * windowsize * 60; %calculate window size in frames
windowaveragetimes = t(windowframes+1:(frame_number-(windowframes+1)));

rp_b = zeros(frame_number, 1);   rp_number = length(rp_locs);   % rp_b will be the binary of the raster rp peaks 

for p = 1:rp_number                                % now its binary
    rp_b(rp_locs(p)) = 1;
end



for i = 1:length(windowaveragetimes)
    windowaverage(i) = mean(rp_b(i:i+windowframes)) - mean(rp_b(i+windowframes:i+2*windowframes));
end

smooth(windowaverage, 10);     % smoothing can change parameter

%% plot result


plot(windowaveragetimes, windowaverage)
ylabel('Moving Window Frequency')
xlabel('Time (Min)'); 

%% save results
save([spath fname '.mat']);    
  
    
%% print result for RP1 081018

cbResult = cb_locs'

rpResult = rp_locs'

end

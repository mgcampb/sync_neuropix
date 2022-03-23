%% SCRIPT FOR SYNCING NEUROPIX SPIKING DATA TO NIDAQ BEHAVIORAL DATA 
% MGC 4/10/2020
% Adapted from Giocomo lab code

% adapted to 2-laser optotagging experiments
% MGC 3/13/2021

%% options
tic
opt = struct;
opt.data_dir = 'D:\neuropix_raw\MC14_20210621_g0'; % dataset to process
opt.dt = 0.001; % temporal bin size for downsampled velocity
opt.lick_thresh = 0.3; % threshold for detecting licks
opt.save_dir = 'D:\neuropix_processed'; % where to save processed data
opt.save_name = 'MC15_20210408.mat'; % file to save processed data to

% experimental parameters (useful to save here)
exp_params = struct;
exp_params.laser1_power = 10; % in mW
exp_params.laser2_power = 10;
exp_params.laser1_target = 'DMS';
exp_params.laser2_target = 'VS';

%% libraries
% need these for loading KS output
addpath(genpath('C:\code\spikes'));
addpath(genpath('C:\code\npy-matlab'));

% need these for processing lick channel
addpath(genpath('C:\code\HGRK_analysis_tools'));

%% location of data
[~,main_name]=fileparts(opt.data_dir);
NIDAQ_file = fullfile(opt.data_dir,strcat(main_name,'_t0.nidq.bin'));
NIDAQ_config = fullfile(opt.data_dir,strcat(main_name,'_t0.nidq.meta'));
spike_dir = fullfile(opt.data_dir,strcat(main_name,'_imec0'));

%% load spike times
fprintf('loading spike data...\n');
sp = loadKSdir(spike_dir);

%% hack for when we entered 30000 into KS2 for NP samp rate instead of true calibrated value (usually something like 30000.27)

% find true NP ap sample rate
np_ap_meta_file = fullfile(spike_dir,strcat(main_name,'_t0.imec0.ap.meta'));
fp_np = fopen(np_ap_meta_file);
dat=textscan(fp_np,'%s %s','Delimiter','=');
names=dat{1};
vals=dat{2};
loc=contains(names,'imSampRate');
true_sampling_rate=str2double(vals{loc});
fclose(fp_np);

% find sample rate that was entered into KS2 (probably 30KHz)
ks_params_file = fullfile(spike_dir,'params.py');
fp_ks = fopen(ks_params_file);
dat=textscan(fp_ks,'%s %s','Delimiter','=');
names=dat{1};
vals=dat{2};
loc=contains(names,'sample_rate');
ks_sampling_rate = str2double(vals{loc});
fclose(fp_ks);

% correct spike times
st=sp.st;
in_samples=st*ks_sampling_rate;
sp.st=in_samples/true_sampling_rate;

%% load nidaq data
fprintf('loading nidaq data...\n');

% get nidaq params
dat=textscan(fopen(NIDAQ_config),'%s %s','Delimiter','=');
names=dat{1};
vals=dat{2};
loc=contains(names,'niSampRate');
daq_sampling_rate=str2double(vals{loc});

nSavedChans = str2double(vals(strcmp(names,'nSavedChans')));

% load nidaq data
fpNIDAQ=fopen(NIDAQ_file);
daq_data=fread(fpNIDAQ,[nSavedChans,Inf],'*int16');
daq_time = (1:size(daq_data,2))/daq_sampling_rate;
fclose(fpNIDAQ);



%% process nidaq data

fprintf('processing nidaq data...\n');
ch_vel = 2;
ch_licks = 3;
ch_rews = 4; % reward valve opening
ch_laser1 = 5; % laser1 trigger
ch_laser2 = 6; % laser2 trigger
ch_cam_sync = 7; % camera sync pulses


% VELOCITY 
% Extract velocity
daq_vel = double(daq_data(ch_vel,:)); % voltage used to compute velocity
daq_vel(daq_vel>prctile(daq_vel,99.9)) = NaN; % remove outliers
daq_vel = fillmissing(daq_vel,'previous'); 
daq_vel = daq_vel-mode(daq_vel);

% downsample velocity
velt = min(daq_time):opt.dt:max(daq_time);
vel = interp1(daq_time,daq_vel,velt,'spline');

% LICKS   
% HyungGoo's lick detection code (modified by Malcolm):
daq_licks = double(daq_data(ch_licks,:));
daq_licks = (daq_licks-min(daq_licks))/(max(daq_licks)-min(daq_licks)); % scale to be 0 to 1, then...
daq_licks = 1-daq_licks; % invert the signal
lick_ts = detect_small_lick_by_deflection_malcolm(daq_licks,opt.lick_thresh,daq_sampling_rate); 
lick_ts = lick_ts/daq_sampling_rate; % convert to seconds

% REWARD
rew_voltage = double(daq_data(ch_rews,:));
% binarize signal
rew_voltage(rew_voltage <= max(rew_voltage)*0.8) = 0;
rew_voltage(rew_voltage > max(rew_voltage)*0.8) = 1;
rew_idx = find(diff(rew_voltage)>0.5)+1;
rew_off_idx = find(diff(rew_voltage)<-0.5)+1;
rew_ts = daq_time(rew_idx)';
rew_off_ts = daq_time(rew_off_idx)';
rew_dur = rew_off_ts-rew_ts;
rew_dur_ms = round(rew_dur*200)*5;

% LASER 1
laser1_voltage = double(daq_data(ch_laser1,:));
% binarize signal
laser1_voltage(laser1_voltage <= max(laser1_voltage)*0.8) = 0;
laser1_voltage(laser1_voltage > max(laser1_voltage)*0.8) = 1;
laser1_onset_idx = find(diff(laser1_voltage)>0.5)+1;
laser1_onset_ts = daq_time(laser1_onset_idx)';
laser1_offset_idx = find(diff(laser1_voltage)<-0.5)+1;
laser1_offset_ts = daq_time(laser1_offset_idx)';
laser1_dur = laser1_offset_ts-laser1_onset_ts;
laser1_dur_ms = round(laser1_dur*1000);

% LASER 2
laser2_voltage = double(daq_data(ch_laser2,:));
% binarize signal
laser2_voltage(laser2_voltage <= max(laser2_voltage)*0.8) = 0;
laser2_voltage(laser2_voltage > max(laser2_voltage)*0.8) = 1;
laser2_onset_idx = find(diff(laser2_voltage)>0.5)+1;
laser2_onset_ts = daq_time(laser2_onset_idx)';
laser2_offset_idx = find(diff(laser2_voltage)<-0.5)+1;
laser2_offset_ts = daq_time(laser2_offset_idx)';
laser2_dur = laser2_offset_ts-laser2_onset_ts;
laser2_dur_ms = round(laser2_dur*1000);

% CAMERA SYNC PULSES
cam_sync_voltage = double(daq_data(ch_cam_sync,:));
% binarize signal
cam_sync_voltage(cam_sync_voltage <= max(cam_sync_voltage)*0.8) = 0;
cam_sync_voltage(cam_sync_voltage > max(cam_sync_voltage)*0.8) = 1;
cam_sync_idx = find(diff(cam_sync_voltage)~=0)+1;
cam_sync_ts = daq_time(cam_sync_idx);

%% CORRECT FOR DRIFT BETWEEN IMEC AND NIDAQ BOARDS

% TWO-PART CORRECTION
% 1. Get sync pulse times relative to NIDAQ and Imec boards.  
% 2. Quantify difference between the two sync pulse times and correct in
% spike.st. 
% PART 1: GET SYNC TIMES RELATIVE TO EACH BOARD
% We already loaded most of the NIDAQ data above. Here, we access the sync
% pulses used to sync Imec and NIDAQ boards together. The times a pulse is
% emitted and registered by the NIDAQ board are stored in syncDatNIDAQ below.

fprintf('correcting drift...\n');
syncDatNIDAQ=daq_data(1,:)>1000;

% convert NIDAQ sync data into time data by dividing by the sampling rate
ts_NIDAQ = strfind(syncDatNIDAQ,[0 1])/daq_sampling_rate; 

% ts_NIDAQ: these are the sync pulse times relative to the NIDAQ board
% Now, we do the same, but from the perspective of the Imec board. 
LFP_config = dir(fullfile(spike_dir,'*.lf.meta'));
LFP_config = fullfile(LFP_config.folder,LFP_config.name);
LFP_file = dir(fullfile(spike_dir,'*.lf.bin'));
LFP_file = fullfile(LFP_file.folder,LFP_file.name);
dat=textscan(fopen(LFP_config),'%s %s','Delimiter','=');
names=dat{1};
vals=dat{2};
loc=contains(names,'imSampRate');
lfp_sampling_rate=str2double(vals{loc});

% for loading only a portion of the LFP data
fpLFP = fopen(LFP_file);
fseek(fpLFP, 0, 'eof'); % go to end of file
fpLFP_size = ftell(fpLFP); % report size of file
fpLFP_size = fpLFP_size/(2*384); 
fclose(fpLFP);

% get the sync pulse times relative to the Imec board
fpLFP=fopen(LFP_file);
fseek(fpLFP,384*2,0);
ftell(fpLFP);
datLFP=fread(fpLFP,[1,round(fpLFP_size/4)],'*int16',384*2); % this step used to take forever
fclose(fpLFP);
syncDatLFP=datLFP(1,:)>10; 
ts_LFP = strfind(syncDatLFP,[0 1])/lfp_sampling_rate;

% ts_LFP: these are the sync pulse times relative to the Imec board
% PART 2: TIME CORRECTION
lfpNIDAQdif = ts_LFP - ts_NIDAQ(1:size(ts_LFP, 2)); % calculate the difference between the sync pulse times
fit = polyfit(ts_LFP, lfpNIDAQdif, 1); % linear fit 
correction_slope = fit(1); % this is the amount of drift we get per pulse (that is, per second)

% save the old, uncorrected data as sp.st_uncorrected and save the new,
% corrected data as sp.st (as many of your analyses are using sp.st).
sp.st_uncorrected = sp.st; % save uncorrected spike times (st)
st_corrected = sp.st - sp.st * correction_slope; % in two steps to avoid confusion
sp.st = st_corrected; % overwrite the old sp.st

%% save processed data
fprintf('saving processed data...\n');
if exist(opt.save_dir,'dir')~=7
    mkdir(opt.save_dir);
end
save(fullfile(opt.save_dir,opt.save_name),...
    'sp','rew_ts','rew_dur_ms','laser1_onset_ts','laser1_offset_ts','laser1_dur_ms',...
    'laser2_onset_ts','laser2_offset_ts','laser2_dur_ms','lick_ts','exp_params',...
    'vel','velt','cam_sync_ts');
fprintf('done.\n');
toc
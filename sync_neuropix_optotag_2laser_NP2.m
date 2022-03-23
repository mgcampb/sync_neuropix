%% SCRIPT FOR SYNCING NEUROPIX SPIKING DATA TO NIDAQ BEHAVIORAL DATA 
% MGC 4/10/2020
% Adapted from Giocomo lab code

% adapted to 2-laser optotagging experiments
% MGC 3/13/2021

% adapted to NP2.0 probes
% sorted using ecephys pipeline
% MGC 6/22/2021

%% options
tic
opt = struct;
opt.run_name = 'MC16_20210626'; % dataset to process
opt.dt = 0.001; % temporal bin size for downsampled velocity
opt.lick_thresh = 0.3; % threshold for detecting licks

% experimental parameters (useful to save here)
exp_params = struct;
exp_params.laser1_power = 10; % in mW
exp_params.laser2_power = 10;
exp_params.laser1_target = 'VS';
exp_params.laser2_target = 'DLS';

% paths
opt.save_dir = 'D:\neuropix_processed'; % where to save processed data
opt.save_name = opt.run_name; % file to save processed data to
opt.ecephys_dir = ['D:\ecephys\catGT\' opt.run_name '\catgt_' opt.run_name '_g0\'];
opt.ks_dir = [opt.ecephys_dir opt.run_name '_g0_imec0\imec0_ks2'];
opt.ni_dir = ['D:\neuropix_raw\' opt.run_name '_g0']; % location of raw nidaq data
opt.sync_ni = [opt.ecephys_dir opt.run_name '_g0_tcat.nidq.XA_0_500.txt'];
opt.sync_np = [opt.ecephys_dir opt.run_name '_g0_imec0\' opt.run_name '_g0_tcat.imec0.ap.SY_384_6_500.txt'];

%% libraries
% need these for loading KS output
addpath(genpath('C:\code\spikes'));
addpath(genpath('C:\code\npy-matlab'));

% need these for processing lick channel
addpath(genpath('C:\code\HGRK_analysis_tools'));

%% location of data
NIDAQ_file = fullfile(opt.ni_dir,strcat(opt.run_name,'_g0_t0.nidq.bin'));
NIDAQ_config = fullfile(opt.ni_dir,strcat(opt.run_name,'_g0_t0.nidq.meta'));

%% load spike times
fprintf('loading spike data...\n');
sp = loadKSdir(opt.ks_dir);

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
% 1. Get sync pulse times from output of catGT
% 2. Map spike times onto NIDAQ stream

fprintf('correcting drift...\n');

% NIDAQ sync pulses:
ts_ni = textscan(fopen(opt.sync_ni),'%f');
ts_np = textscan(fopen(opt.sync_np),'%f');
ts_ni = ts_ni{1};
ts_np = ts_np{1};

assert(numel(ts_ni)==numel(ts_np));

% PART 2: TIME CORRECTION
ts_diff = ts_np - ts_ni; % calculate the difference between the sync pulse times
fit = polyfit(ts_np, ts_diff, 1); % linear fit 
correction_slope = fit(1); % this is the amount of drift we get per pulse (that is, per second)

% save the old, uncorrected data as sp.st_uncorrected and save the new,
% corrected data as sp.st (as many of your analyses are using sp.st).
sp.st_uncorrected = sp.st; % save uncorrected spike times (st)
st_corrected = sp.st - sp.st * correction_slope; % in two steps to avoid confusion
sp.st = st_corrected; % overwrite the old sp.st

%% load additional information about units: QC metrics, waveform metrics, channel map
sp.chan_map = readNPY(fullfile(opt.ks_dir,'channel_map.npy'));
sp.waveform_metrics = readtable(fullfile(opt.ks_dir,'waveform_metrics.csv'));
sp.waveform_metrics = sp.waveform_metrics(:,2:end);
sp.qc_metrics = readtable(fullfile(opt.ks_dir,'metrics.csv'));
sp.qc_metrics = sp.qc_metrics(:,2:end);
sp.mean_waveforms = readNPY(fullfile(opt.ks_dir,'mean_waveforms.npy'));

% select cells that were manually labeled
idx = nan(size(sp.cgs));
for i = 1:numel(idx)
    idx_this = find(sp.qc_metrics.cluster_id==sp.cids(i));
    if ~isempty(idx_this)
        idx(i) = idx_this;
    end
end
idx = idx(~isnan(idx));
sp.waveform_metrics = sp.waveform_metrics(idx,:);
sp.qc_metrics = sp.qc_metrics(idx,:);
sp.mean_waveforms = sp.mean_waveforms(sp.cids(ismember(sp.cids,sp.qc_metrics.cluster_id))+1,:,:);
sp.mean_waveforms_cluster_id = sp.cids(ismember(sp.cids,sp.qc_metrics.cluster_id));

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
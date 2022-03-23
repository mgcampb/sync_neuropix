%% SCRIPT FOR SYNCING NEUROPIX SPIKING DATA TO BPOD BEHAVIORAL DATA
% MGC 9/2/2021

% adapted to concatenate two Bpod protocols (for Sara's experiments)
% 9/27/2021

%% options
tic
opt = struct;
opt.run_name = 'SM138_20211005'; % dataset to process

strspl = strsplit(opt.run_name,'_');
mouse = strspl{1};
date = strspl{2};

% experimental parameters (useful to save here)
exp_params = struct;
exp_params.laser_power = [10]; % in mW
exp_params.laser_target = {'VS'};

% paths
opt.save_dir = 'F:\neuropix_processed'; % where to save processed data
opt.save_name = opt.run_name; % file to save processed data to
opt.ecephys_dir = ['F:\ecephys\catGT\catgt_' opt.run_name '_g0\'];
opt.ks_dir = [opt.ecephys_dir opt.run_name '_g0_imec0\imec0_ks3_thresh99'];
opt.bpod_dir1 = ['F:\Bpod_Data\' mouse '\Eshel2016\Session Data']; % location of raw bpod data (1/2)
opt.bpod_dir2 = ['F:\Bpod_Data\' mouse '\OptoTagging\Session Data']; % location of raw bpod data (2/2)
opt.sync_imec = [opt.ecephys_dir opt.run_name '_g0_imec0\' opt.run_name '_g0_tcat.imec0.ap.SY_384_6_0.txt'];

%% libraries
% need these for loading KS output
addpath(genpath('C:\code\spikes'));
addpath(genpath('C:\code\npy-matlab'));

% need these for processing lick channel
addpath(genpath('C:\code\HGRK_analysis_tools'));

%% load spike times
fprintf('loading spike data...\n');
sp = loadKSdir(opt.ks_dir);

%% load sync data (Neuropix)
ts_imec = textscan(fopen(opt.sync_imec),'%f');
ts_imec = ts_imec{1};

%% load Bpod data from both sessions

fn_bpod1 = dir(opt.bpod_dir1);
fn_bpod1 = {fn_bpod1.name}';
fn_bpod1 = fn_bpod1{contains(fn_bpod1,date)};
bpod_dat1 = load(fullfile(opt.bpod_dir1,fn_bpod1));

fn_bpod2 = dir(opt.bpod_dir2);
fn_bpod2 = {fn_bpod2.name}';
fn_bpod2 = fn_bpod2{contains(fn_bpod2,date)};
bpod_dat2 = load(fullfile(opt.bpod_dir2,fn_bpod2));

num_trial_protocol1 = numel(bpod_dat1.SessionData.TrialStartTimestamp);

ts_imec = ts_imec([1:num_trial_protocol1 num_trial_protocol1+2:end]);

offset = ts_imec(num_trial_protocol1+1)-ts_imec(1);
bpod_dat2.SessionData.TrialStartTimestamp = bpod_dat2.SessionData.TrialStartTimestamp + offset;
bpod_dat2.SessionData.TrialEndTimestamp = bpod_dat2.SessionData.TrialEndTimestamp + offset;

ts_bpod = [bpod_dat1.SessionData.TrialStartTimestamp'+0.05; ...
    bpod_dat2.SessionData.TrialStartTimestamp'];
assert(numel(ts_bpod)==numel(ts_imec));

SessionData1 = bpod_dat1.SessionData;
SessionData2 = bpod_dat2.SessionData;

%% Plot sync pulse difference
% should look like a line with one outlier corresponding to the gap between
% BpodSession1 and BpodSession2
figure; hold on
title([mouse ' ' date])
plot([0 max(diff(ts_imec))],[0 max(diff(ts_imec))])
plot(diff(ts_imec),diff(ts_bpod),'ko');
xlabel('sync pulse difference, IMEC (sec)');
ylabel('sync pulse difference, Bpod (sec)');

%% Sync spike times to Bpod
sp.st_uncorrected = sp.st;
sp.st = interp1(ts_imec,ts_bpod,sp.st);

%% Remove unnecessary data
sp = rmfield(sp,'spikeTemplates');
sp = rmfield(sp,'tempScalingAmps');
sp = rmfield(sp,'temps');
sp = rmfield(sp,'winv');
sp = rmfield(sp,'pcFeat');
sp = rmfield(sp,'pcFeatInd');

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
    'sp','SessionData1','SessionData2','exp_params');
fprintf('done.\n');
toc
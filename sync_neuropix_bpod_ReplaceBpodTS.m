%% SCRIPT FOR SYNCING NEUROPIX SPIKING DATA TO BPOD BEHAVIORAL DATA
% MGC 9/2/2021

% add NIDAQ data (e.g. running speed)
% script will detect if it's there and skip if it's not
% Added linfit/quadfit synchronization to get whole session
% MGC 3/23/2022

% New sync approach:
% Just replace Bpod trial start time stamps with Imec's
% (Simpler, and we keep spikes for the whole session)
% (Drift should be minimal within a trial)

%% options
tic
opt = struct;
% datasets to process:
opt.save_dir_orig = 'I:\My Drive\UchidaLab\DA_independence\neuropix_processed\ks3_thresh99\';

opt.run_name = dir(fullfile(opt.save_dir_orig,'*.mat'));
opt.run_name = {opt.run_name.name}';
for i = 1:numel(opt.run_name)
    opt.run_name{i} = opt.run_name{i}(1:end-4);
end

keep = ~strcmp(opt.run_name,'MC34_20220110');
opt.run_name = opt.run_name(keep);

opt.nidaq_tbin = 0.001; % in seconds
opt.nidaq_ch = [2]; % which nidaq channels to keep
opt.nidaq_name = {'velocity'}; % names for the above channels

% % experimental parameters (useful to save here)
% exp_params = struct;
% exp_params.bpod_protocol = 'OdorLaser';
% exp_params.laser_power = [10,10,10]; % in mW
% exp_params.laser_target = {'VS','DMS','DLS'};

%% libraries
% need these for loading KS output
addpath(genpath('C:\code\spikes'));
addpath(genpath('C:\code\npy-matlab'));

% need these for loading meta data
addpath(genpath('C:\code\SGLX_tools'));

% need these for processing lick channel
addpath(genpath('C:\code\HGRK_analysis_tools'));

for sesh_idx = 1:numel(opt.run_name)
    
    run_name = opt.run_name{sesh_idx};
    strspl = strsplit(run_name,'_');
    mouse = strspl{1};
    date = strspl{2};

    dat_orig = load(fullfile(opt.save_dir_orig,run_name));
    exp_params = dat_orig.exp_params;

    % paths
    opt.save_dir = 'I:\My Drive\UchidaLab\DA_independence\neuropix_processed\replace_bpod_ts'; % where to save processed data
    opt.save_name = run_name; % file to save processed data to
    opt.ecephys_dir = ['F:\ecephys\catGT\catgt_' run_name '_g0\'];
    opt.ks_dir = [opt.ecephys_dir run_name '_g0_imec0\imec0_ks3'];
    opt.bpod_dir = ['F:\Bpod_Data\' mouse '\' exp_params.bpod_protocol '\Session Data']; % location of raw bpod data
    opt.sync_imec = [opt.ecephys_dir run_name '_g0_imec0\' run_name '_g0_tcat.imec0.ap.SY_384_6_0.txt'];
    opt.ni_dir = ['F:\neuropix_raw\' run_name '_g0\'];

    %% load spike times
    fprintf('loading spike data...\n');
    sp = loadKSdir(opt.ks_dir);

    %% load sync data (Neuropix)
    ts_imec = textscan(fopen(opt.sync_imec),'%f');
    ts_imec = ts_imec{1};

    %% load Bpod data
    fn_bpod = dir(opt.bpod_dir);
    fn_bpod = {fn_bpod.name}';
    fn_bpod = fn_bpod{contains(fn_bpod,date)};
    load(fullfile(opt.bpod_dir,fn_bpod));

    ts_bpod = SessionData.TrialStartTimestamp';
    if numel(ts_bpod)~=numel(ts_imec) % sometimes missing one sync pulse (early mistakes)
        if strcmp(run_name,'MC34_20220111') % weird session where forgot to connect sync cable
            ts_imec = ts_imec(2:end);
        else
            ts_imec = ts_imec(1:399);
        end
    end
    assert(corr(diff(ts_bpod),diff(ts_imec))>0.95);

    % New sync approach: Just replace Bpod trial start timestamps with Imec's
    % (Drift should be minimal within a trial)
    trial_diff = SessionData.TrialEndTimestamp - SessionData.TrialStartTimestamp;
    SessionData.TrialStartTimestamp = ts_imec';
    SessionData.TrialEndTimestamp = ts_imec' + trial_diff;

    %% Old sync approach: Sync spike times to Bpod
    % sp.st_uncorrected = sp.st; % save original spike times

    % % interp1. Add 10 seconds to last ts to get the last trial.   
    % sp.st = interp1([ts_imec; ts_imec(end)+10],[ts_bpod; ts_bpod(end)+10],sp.st_uncorrected); 

    %% load NIDAQ data
    NidaqData = [];
    ni_bin = dir(fullfile(opt.ni_dir,'*.bin'));

    if ~isempty(ni_bin)
        ni_bin = ni_bin.name;
        ni_meta = ReadMeta(ni_bin,opt.ni_dir);
        ni_nSavedChans = str2double(ni_meta.nSavedChans);
        ni_SampRate = str2double(ni_meta.niSampRate);
    
        fprintf('loading nidaq data...\n');
        fpNIDAQ = fopen(fullfile(opt.ni_dir,ni_bin));
        daq_data = fread(fpNIDAQ,[ni_nSavedChans,Inf],'*int16');
        fclose(fpNIDAQ);
    
        % find rising edges of sync pulses on first NI channel
        daq_time = (0:size(daq_data,2)-1)/ni_SampRate;
        ni_sync = daq_data(1,:);
        ni_sync = (ni_sync-min(ni_sync))/(max(ni_sync)-min(ni_sync)); % normalize between 0 and 1
        ni_sync = round(ni_sync); % get rid of any in between points
        ni_sync_idx = find(diff(ni_sync)==1)+1;
        ts_ni = daq_time(ni_sync_idx)';
    
        if numel(ts_ni)~=numel(ts_imec) % sometimes missing one sync pulse (early mistakes)
            if strcmp(run_name,'MC34_20220111') % weird session where forgot to connect sync cable
                ts_ni = ts_ni(2:end);
            else
                ts_ni = ts_ni(1:399);
            end
        end
        assert(corr(diff(ts_ni),diff(ts_imec))>0.95);

        % sync to Imec
        fit = polyfit(ts_ni,ts_imec,1); % linear fit
        daq_time_transf = daq_time * fit(1) + fit(2);
        daqt = min(daq_time_transf):opt.nidaq_tbin:max(daq_time_transf);

        NidaqData = struct;
        NidaqData.daqt = daqt;
        NidaqData.data = nan(numel(opt.nidaq_ch),numel(daqt));
        for daq_ch = 1:numel(opt.nidaq_ch)
            NidaqData.data(daq_ch,:) = interp1(daq_time_transf,double(daq_data(opt.nidaq_ch(daq_ch),:)),daqt);
        end
        NidaqData.name = opt.nidaq_name;
    end 

    

    %% Remove unnecessary data
    sp = rmfield(sp,'spikeTemplates');
    sp = rmfield(sp,'tempScalingAmps');
    sp = rmfield(sp,'temps');
    sp = rmfield(sp,'winv');
    sp = rmfield(sp,'pcFeat');
    sp = rmfield(sp,'pcFeatInd');

    %% Mean waveforms and metrics from original files
    sp.mean_waveforms = dat_orig.sp.mean_waveforms;
    sp.metrics = dat_orig.sp.metrics;

    %% save processed data
    fprintf('saving processed data...\n');
    if exist(opt.save_dir,'dir')~=7
        mkdir(opt.save_dir);
    end
    save(fullfile(opt.save_dir,opt.save_name),...
        'sp','SessionData','exp_params','NidaqData');
    fprintf('done.\n');
    toc
end
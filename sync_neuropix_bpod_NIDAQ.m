%% SCRIPT FOR SYNCING NEUROPIX SPIKING DATA TO BPOD BEHAVIORAL DATA
% MGC 9/2/2021

% add NIDAQ data (e.g. running speed)
% MGC 3/23/2022

%% options
tic
opt = struct;
% datasets to process:
opt.save_dir_orig = 'D:\neuropix_processed';

opt.run_name = dir(fullfile(opt.save_dir_orig,'*.mat'));
opt.run_name = {opt.run_name.name}';
for i = 1:numel(opt.run_name)
    opt.run_name{i} = opt.run_name{i}(1:end-4);
end
opt.run_name = opt.run_name(contains(opt.run_name,'MC25'));

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
    opt.save_dir = 'C:\data\neuropix_processed_quadfit'; % where to save processed data
    opt.save_name = run_name; % file to save processed data to
    opt.ecephys_dir = ['D:\ecephys\catGT\catgt_' run_name '_g0\'];
    opt.ks_dir = [opt.ecephys_dir run_name '_g0_imec0\imec0_ks3'];
    opt.bpod_dir = ['D:\Bpod_Data\' mouse '\' exp_params.bpod_protocol '\Session Data']; % location of raw bpod data
    opt.sync_imec = [opt.ecephys_dir run_name '_g0_imec0\' run_name '_g0_tcat.imec0.ap.SY_384_6_0.txt'];
    opt.ni_dir = ['D:\neuropix_raw\' run_name '_g0\'];

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
    assert(numel(ts_bpod)==numel(ts_imec));
    assert(corr(diff(ts_bpod),diff(ts_imec))>0.95);
    
    %% load NIDAQ data
    ni_bin = dir(fullfile(opt.ni_dir,'*.bin'));
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

    %% Sync spike times to Bpod
    sp.st_uncorrected = sp.st; % save original spike times

    % Old way: interp1. Add 10 seconds to last ts to get the last trial.   
    sp.st_interp = interp1([ts_imec; ts_imec(end)+10],[ts_bpod; ts_bpod(end)+10],sp.st_uncorrected); 

    % New way: quadratic fit. Gets the whole session, including parts with no sync pulses.
    fit = polyfit(ts_imec,ts_bpod,2); 
    sp.st = sp.st_uncorrected.^2 * fit(1) + sp.st_uncorrected * fit(2) + fit(3);

    %% Sync NIDAQ data to Bpod: Velocity

    % get velocity
    velt_orig = 0:0.001:max(daq_time); % 1ms time bins
    vel_orig = interp1(daq_time,double(daq_data(2,:)),velt_orig);  

    fit = polyfit(ts_ni,ts_bpod,2); % quadratic fit
    velt_transf = velt_orig.^2 * fit(1) + velt_orig * fit(2) + fit(3);
    velt = min(velt_transf):0.001:max(velt_transf);
    vel = interp1(velt_transf,vel_orig,velt);
    
    % lightly smooth velocity
    win = gausswin(50); % 50 ms gaussian window (not sure of exact SD)
    win = win/sum(win);
    vel = conv(vel,win,"same");

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

    %% load additional information about units: QC metrics, waveform metrics, channel map
    
    % sp.chan_map = readNPY(fullfile(opt.ks_dir,'channel_map.npy'));
    % sp.waveform_metrics = readtable(fullfile(opt.ks_dir,'waveform_metrics.csv'));
    % sp.waveform_metrics = sp.waveform_metrics(:,2:end);
    % sp.qc_metrics = readtable(fullfile(opt.ks_dir,'metrics.csv'));
    % sp.qc_metrics = sp.qc_metrics(:,2:end);
    % sp.mean_waveforms = readNPY(fullfile(opt.ks_dir,'mean_waveforms.npy'));
    % 
    % % select cells that were manually labeled
    % idx = nan(size(sp.cgs));
    % for i = 1:numel(idx)
    %     idx_this = find(sp.qc_metrics.cluster_id==sp.cids(i));
    %     if ~isempty(idx_this)
    %         idx(i) = idx_this;
    %     end
    % end
    % idx = idx(~isnan(idx));
    % sp.waveform_metrics = sp.waveform_metrics(idx,:);
    % sp.qc_metrics = sp.qc_metrics(idx,:);
    % sp.mean_waveforms = sp.mean_waveforms(sp.cids(ismember(sp.cids,sp.qc_metrics.cluster_id))+1,:,:);
    % sp.mean_waveforms_cluster_id = sp.cids(ismember(sp.cids,sp.qc_metrics.cluster_id));

    %% remove double-counted spikes
    % Not necessary - already none - MGC 11/9/2021

    % fprintf('removing double-counted spikes...\n');
    % good_cells = sp.cids(sp.cgs==2);
    % st_orig = sp.st(ismember(sp.clu,good_cells));
    % clu_orig = sp.clu(ismember(sp.clu,good_cells));
    % 
    % st_uniq = [];
    % clu_uniq = [];
    % for i = 1:numel(good_cells)
    %     st_this = unique(st_orig(clu_orig==good_cells(i)));
    %     st_uniq = [st_uniq; st_this];
    %     clu_uniq = [clu_uniq; good_cells(i)*ones(numel(st_this),1)];
    % end
    % 
    % [st_uniq,sort_idx] = sort(st_uniq);
    % clu_uniq = clu_uniq(sort_idx);
    % sp.st_uniq = st_uniq;
    % sp.clu_uniq = clu_uniq;
    % 
    % fprintf('\t removed %d double-counted spikes\n',numel(st_orig)-numel(st_uniq));

    %% save processed data
    fprintf('saving processed data...\n');
    if exist(opt.save_dir,'dir')~=7
        mkdir(opt.save_dir);
    end
    save(fullfile(opt.save_dir,opt.save_name),...
        'sp','SessionData','exp_params','velt','vel');
    fprintf('done.\n');
    toc
end
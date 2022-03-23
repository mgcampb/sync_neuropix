%% SCRIPT FOR SYNCING NEUROPIX SPIKING DATA TO BPOD BEHAVIORAL DATA
% MGC 9/2/2021

% add NIDAQ data (e.g. running speed)
% MGC 3/23/2022

%% options
tic
opt = struct;
% datasets to process:
opt.run_name = {'MC34_20220111'};

% experimental parameters (useful to save here)
exp_params = struct;
exp_params.bpod_protocol = 'OdorLaser';
exp_params.laser_power = [10,10,10]; % in mW
exp_params.laser_target = {'VS','DMS','DLS'};

%% libraries
% need these for loading KS output
addpath(genpath('C:\code\spikes'));
addpath(genpath('C:\code\npy-matlab'));

% need these for processing lick channel
addpath(genpath('C:\code\HGRK_analysis_tools'));

for sesh_idx = 1:numel(opt.run_name)
    
    run_name = opt.run_name{sesh_idx};
    strspl = strsplit(run_name,'_');
    mouse = strspl{1};
    date = strspl{2};

    % paths
    opt.save_dir = 'F:\neuropix_processed'; % where to save processed data
    opt.save_name = run_name; % file to save processed data to
    opt.ecephys_dir = ['F:\ecephys\catGT\catgt_' run_name '_g0\'];
    opt.ks_dir = [opt.ecephys_dir run_name '_g0_imec0\imec0_ks3'];
    opt.bpod_dir = ['F:\Bpod_Data\' mouse '\' exp_params.bpod_protocol '\Session Data']; % location of raw bpod data
    opt.sync_imec = [opt.ecephys_dir run_name '_g0_imec0\' run_name '_g0_tcat.imec0.ap.SY_384_6_0.txt'];
    opt.nidaq = ['F:\neuropix_raw\' run_name '_g0\'];

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
        'sp','SessionData','exp_params');
    fprintf('done.\n');
    toc
end
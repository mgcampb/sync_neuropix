%% SCRIPT FOR SYNCING NEUROPIX SPIKING DATA TO NIDAQ BEHAVIORAL DATA 
% MGC 4/10/2020
% Adapted from Giocomo lab code

% Adapted for OptoInhibition test, MGC 1/6/2022

%% options
tic
opt = struct;
opt.session = 'MC45_20220222_OptoInhibitionTest_Jaws'; % dataset to process
opt.main_name = [opt.session '_g0'];
opt.data_dir_spike = ['F:\ecephys\catGT\catgt_' opt.main_name '\' opt.main_name '_imec0\']; % location of catGT-processed raw data
opt.data_dir_ks = [opt.data_dir_spike '\imec0_ks3']; % location of KS output
opt.data_dir_nidaq = ['F:\neuropix_raw\' opt.main_name '\']; % location of nidaq data
opt.save_dir = 'F:\neuropix_processed'; % where to save processed data
opt.save_name = opt.session; % file to save processed data to

%% libraries
% need these for loading KS output
addpath(genpath('C:\code\spikes'));
addpath(genpath('C:\code\npy-matlab'));

% need these for processing lick channel
addpath(genpath('C:\code\HGRK_analysis_tools'));

%% location of NIDAQ data
NIDAQ_file = fullfile(opt.data_dir_nidaq,strcat(opt.main_name,'_t0.nidq.bin'));
NIDAQ_config = fullfile(opt.data_dir_nidaq,strcat(opt.main_name,'_t0.nidq.meta'));

%% load spike times
fprintf('loading spike data...\n');
sp = loadKSdir(opt.data_dir_ks);

%% hack for when we entered 30000 into KS2 for NP samp rate instead of true calibrated value (usually something like 30000.27)

% find true NP ap sample rate
np_ap_meta_file = fullfile(opt.data_dir_spike,strcat(opt.main_name,'_tcat.imec0.ap.meta'));
fp_np = fopen(np_ap_meta_file);
dat=textscan(fp_np,'%s %s','Delimiter','=');
names=dat{1};
vals=dat{2};
loc=contains(names,'imSampRate');
true_sampling_rate=str2double(vals{loc});
fclose(fp_np);

% find sample rate that was entered into KS (probably 30KHz)
ks_params_file = fullfile(opt.data_dir_ks,'params.py');
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
fpNIDAQ=fopen(NIDAQ_file);
daq_data=fread(fpNIDAQ,[9,Inf],'*int16');
fclose(fpNIDAQ);

% get the nidaq sample rate
dat=textscan(fopen(NIDAQ_config),'%s %s','Delimiter','=');
names=dat{1};
vals=dat{2};
loc=contains(names,'niSampRate');
daq_sampling_rate=str2double(vals{loc});


%% process nidaq data

fprintf('processing nidaq data...\n');

ch_trig = 5;
ch_power = 8;

% downsample DAQ data to 1 kHz
daq_time_orig = (0:size(daq_data,2)-1)/daq_sampling_rate;
daq_time = 0:0.001:max(daq_time_orig);

laser_trig = interp1(daq_time_orig,double(daq_data(ch_trig,:)),daq_time);
laser_power = interp1(daq_time_orig,double(daq_data(ch_power,:)),daq_time);

% binarize the trigger
laser_trig(laser_trig<max(laser_trig)/2) = 0;
laser_trig(laser_trig>max(laser_trig)/2) = 1;

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
fpNIDAQ=fopen(NIDAQ_file);
daq_data=fread(fpNIDAQ,[9,Inf],'*int16');
syncDatNIDAQ=daq_data(1,:)>1000;
ts_NIDAQ = strfind(syncDatNIDAQ,[0 1])/daq_sampling_rate; 

% ts_NIDAQ: these are the sync pulse times relative to the NIDAQ board
% Now, we do the same, but from the perspective of the Imec board. 
LFP_config = dir(fullfile(opt.data_dir_spike,'*.lf.meta'));
LFP_config = fullfile(LFP_config.folder,LFP_config.name);
LFP_file = dir(fullfile(opt.data_dir_spike,'*.lf.bin'));
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
    'sp','daq_time','laser_trig','laser_power');
fprintf('done.\n');
toc
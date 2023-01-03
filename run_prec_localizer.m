%% start

% addpath C:\toolbox\fieldtrip-20220707
% ft_defaults


clear; close all; clc;

% set directories
lochome = getuserdir();
if strcmp(lochome(1),'C'),
    lochome = 'C:\';
end
locproj = fullfile(lochome,'vcnl','projects');
loccurp = fullfile(locproj,'eeg-2step-replay');
locrawd = fullfile(loccurp,'raw');
locderi = fullfile(loccurp,'derivatives');
locexpr = fullfile(loccurp,'experiment');
locruns = fullfile(loccurp,'run');
locthis = fullfile(loccurp,'run','preprocessing');

sub_id = 'sub-01';
filepath = fullfile(locrawd,sub_id);

%%


file_eeg = [sub_id '_2step-replay_localizer.eeg'];
file_evt = [sub_id '_events-bst.mat'];
file_pro = [sub_id '_projection-ssp.mat'];

% read continuous data after filtering
% data_eeg.trial has the continuous eeg
cfg = [];
cfg.lpfilter   =  'yes';
cfg.lpfreq     =  40;
cfg.dataset = fullfile(filepath,file_eeg);
data_eeg  = ft_preprocessing(cfg);


% %% option setting
% docorrect = 1;
% 
% %% need to compare the performance of data cleaning 
% 
% % eeg correction based on the SSP (signal-space projection)
% projector = load(fullfile(filepath,file_pro));
% 
% % projection rule
% % W: signal space
% % I: ssp or independent component
% % S: eeg signal
% 
% % W*I = S
% % I = pinv(W)*S;
% % W_ * I = S_ : corrected signal
% 
% %%%%%%%%%%%%%% SET %%%%%%%%%%%%%%%%%
% cur_proj = 1;
% W = projector.Projector(cur_proj).Components;
% iW = pinv(W);
% M = projector.Projector(cur_proj).CompMask;
% W_ = W;
% W_(:,find(M)) = 0;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% data_clean = data_eeg;
% 
% if docorrect
%     temp = data_clean.trial{1};
%     I = iW*temp;
%     temp = W_*I;
%     data_clean.trial = {temp};
% end

%%
% rereferencing with the linked mastoid
% compared to data_eeg, LM data is added into 32nd channel
cfg = [];
cfg.reref = 'yes';
cfg.channel = 'all';
cfg.implicitref = 'LM'; % once implicit channel is delcared
cfg.refchannel = {'LM','RM'};


% data_clean_lm = ft_preprocessing(cfg,data_clean);
data_eeg_lm = ft_preprocessing(cfg,data_eeg);


%% define trials
start_t = -0.2;
end_t = 2.0;

% first choice_1
cfg = [];
cfg.dataset = fullfile(filepath,file_eeg);
cfg.trialfun = 'trialfunc'; % 항상 있어야
cfg.trialdef.prestim = start_t;
cfg.trialdef.poststim = end_t;
cfg.trialdef.eventtype = 'Stimulus';
% cfg.trialdef.eventvalue = {'S  1'; 'S 12'; 'S 21'; 'S107'; 'S108'; 'S112'; 'S121'; 'S140'; 'S141'; 'S142'; 'S203'; 'S204'; 'S205'; 'S206'; 'S213';'S214'; 'S215'; 'S216'; 'S240'; 'S243'; 'S244'; 'S245'; 'S246'; 'S250'; 'S251' }';
cfg.trialdef.eventvalue = {'S 12'; 'S 21'; 'S112'; 'S121'};
trig = ft_definetrial(cfg);
segment = ft_redefinetrial(trig, data_eeg_lm);

% cfg.trials.eventvalue = 'all';
% trig2 = ft_definetrial(cfg);
% segment2 = ft_redefinetrial(trig2, data_eeg_lm);


%% select trials for artifact rejection     
% choose trials meaning the first choice only;
cfg = [];
% cfg.trials = segment.trialinfo == 2 | segment.trialinfo == 3 | segment.trialinfo == 6 | segment.trialinfo == 7;
cfg.trials = segment.trialinfo == 1 | segment.trialinfo == 2 | segment.trialinfo == 3 | segment.trialinfo == 4;
data_f_clean = ft_selectdata(cfg, segment);


%% reject trials
% 0: without artifact / 1: with artifact
    
% load artifact file from brainstorm
arti_mat = load(fullfile(filepath,file_evt)).events;

% blink signal(1~7hz) rejection
% find trials overlapping with blink point
blink_p = arti_mat(27).times;
blink_temp = [blink_p-0.5; blink_p+0.5]';
blink_in = [];
for i=1:length(blink_temp)
    blink_t =[];
    for j=1:46
        blink_t = [blink_t; blink_temp(i,1)+(j-1)*0.1]; %%??
    end
    blink_in = [blink_in; blink_t'];
end

blink = [];
for d=1:length(data_f_clean.sampleinfo)
    searching = 0;
    check = 0;
    while(~searching)
        for b=1:length(blink_p)
            for in=1:46 % the number of elements in one row
                if data_f_clean.sampleinfo(d,1)/500 <= blink_in(b,in) && blink_in(b,in) <= data_f_clean.sampleinfo(d,2)
                    blink = [blink; 1];
                    check = 1;
                    break;
                end
            end
            if check == 1
                break;
            elseif data_f_clean.sampleinfo(d,2)/500 < blink_in(b+1,1)
                blink = [blink; 0];
                break;
            elseif b == length(blink_p)
                blink = [blink; 0];
                break;
            end
        end
        searching = 1;
    end
end

% fast noise signal rejection
% find trials overlapping with fast noise
high_freq_p = arti_mat(28).times;
high_freq_temp = high_freq_p';
high_freq = [];

for d=1:length(data_f_clean.sampleinfo)
    searching = 0;
    check = 0;
    while(~searching)
        for b=1:length(high_freq_temp)
            % make interval vector per intervals
            high_freq_n = [];
            high_freq_n = high_freq_temp(b,1):0.1:high_freq_temp(b,2);
            if b ~= length(high_freq_temp)
                high_freq_post = [];
                high_freq_post = high_freq_temp(b+1,1):0.1:high_freq_temp(b+1,2);
            end
            for in=1:length(high_freq_n) % the number of elements in one row
                if data_f_clean.sampleinfo(d,1)/500 <= high_freq_n(in) && high_freq_n(in) <= data_f_clean.sampleinfo(d,2)
                    high_freq = [high_freq; 1];
                    check = 1;
                    break;
                end
            end
            if check == 1
                break;
            elseif data_f_clean.sampleinfo(d,2)/500 < high_freq_post(1)
                high_freq = [high_freq; 0];
                break;
            elseif b == length(high_freq_temp)
                high_freq = [high_freq; 0];
                break;
            end
        end
        searching = 1;
    end
end

% select clean trials

clean_t = [];
art_t = [blink high_freq];
for i=1:length(art_t)    
    check = sum(art_t(i,:));
    if check == 0
        clean_t = [clean_t; 0];
    else
        clean_t = [clean_t; 1];
    end
end

data_f_clean.clean_trial = clean_t;

% clean data
trials=[];
cfg = [];
cfg.trials = data_f_clean.clean_trial == 0; % 하나 빼고 다 1 인데요..?
clean_data = ft_selectdata(cfg, data_f_clean);
trials = [trials; length(clean_data.trial)];

%% calculate against mean signal between -0.2 and 0  
% frequent event
trig.demean = 'yes';
trig.baselinewindow = [-0.2 0];
    
cut_data = ft_preprocessing(trig, clean_data);
    
%% segment data 
% first choice
cfg = [];
% cfg.trials = segment.trialinfo == 2 | segment.trialinfo == 3 | segment.trialinfo == 6 | segment.trialinfo == 7;
cfg.trials = cut_data.trialinfo == 1 | cut_data.trialinfo == 2 | cut_data.trialinfo == 3 | cut_data.trialinfo == 4;
data_fc = ft_selectdata(cfg, cut_data);

%% run time-lock analysis
cfg = [];
fc_an = ft_timelockanalysis(cfg, data_fc);


%% run grand-mean ERP calculation 
cfg = [];
cfg.method = 'within';
    
grand_fc = ft_timelockgrandaverage(cfg, fc_an);



%% plot results
cfg=[];
cfg.channel = 'Pz'%{'Fz' 'Pz' 'Cz'}; %pick!
% cfg.channel = 'all';
figure;
ERP = ft_singleplotER(cfg, grand_fc);
    
%% save data
filepath = fullfile(locderi,sub_id);
% set save folder path
save_path = fullfile(filepath, [sub_id '_localizer_preprocessing_SSP.mat']);
ERP_path = fullfile(filepath, [sub_id '_localizer_ERP_SSP.jpg']);
    
% save data
save(save_path, 'grand_fc');
saveas(gcf, ERP_path);
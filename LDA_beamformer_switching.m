% clear all
% global ft_default
% ft_default.spmversion = 'spm12'; % Force SPM12, SPM8 doesn't go well with mac + 2017b
% ft_defaults % This loads the rest of the defaults
% 
% files = dir('*.mat');
% 
% for i=1:length(files)
% 
%     load(files(i).name)
% stay_trials   = [events_allBlocks.cuechstay;events_allBlocks.cueenstay];
% switch_trials = [events_allBlocks.cuechswitch;events_allBlocks.cueenswitch];
% 
% ncond_1 = numel(stay_trials);
% ncond_2 = numel(switch_trials);
% 
% cfg         = [];
% cfg.trials  = stay_trials;
% cfg.latency = [-0.1 0.8];
% allstay     = ft_selectdata(cfg,all_blocks_clean);
% 
% cfg         = [];
% cfg.trials  = switch_trials;
% cfg.latency = [-0.1 0.8];
% allswitch   = ft_selectdata(cfg,all_blocks_clean);
% allswitch_ave=ft_timelockanalysis([],allswitch);
% allstay_ave=ft_timelockanalysis([],allstay);
% 
% cfg=[];
% cfg.operation='subtract';
% cfg.parameter='avg';
% stayvsswitch=ft_math(cfg,allstay_ave,allswitch_ave);
% 
% % allswitch_aveall{i} = allswitch_ave;
% % allstay_aveall{i} = allstay_ave;
% 
% ival   = [0.15 .25];
% topo_dat = stayvsswitch;
% time_idx = (topo_dat.time >= ival(1) & topo_dat.time <= ival(2));
% P        = mean(topo_dat.avg(:, time_idx),2);
% %% LDA beamforming ------------------------------------------------
% % To obtain the LDA beamformer, we need the spatial pattern P and the
% % corresponding continuous or epoched data (NOT the averaged data).
% [w,switchlda,C] = LDAbeamformer(P,allswitch_ave,.1);
% [w,staylda,C] = LDAbeamformer(P,allstay_ave,.1);
% 
% %% Plot covariance matrix, spatial pattern P, and spatial filter w
% % Get channel layout for plotting
% % cfg              = [];
% % cfg.grad         = topo_dat.grad;
% % cfg.skipscale    = 'yes';
% % cfg.skipcomnt    = 'yes';
% % cfg.grad.balance = rmfield(cfg.grad.balance,'reject');
% % lay              = ft_prepare_layout(cfg);
% %
% % % Use the low-level plotting function for plotting arbitrary vectors
% % cfg_topo = {'mask' lay.mask, 'datmask' [], 'interplim' 'mask'};
% % cfg_lay  = { 'point' 'no' 'box' 'no' 'label' 'no' };
% %
% % figure
% % subplot(1,3,1),
% % imagesc(C),title('Covariance matrix')
% % subplot(1,3,2),
% % ft_plot_topo(lay.pos(:,1), lay.pos(:,2), P, cfg_topo{:});
% % ft_plot_lay(lay, cfg_lay{:})
% % title('Spatial pattern P')
% % subplot(1,3,3),
% % ft_plot_topo(lay.pos(:,1), lay.pos(:,2), w, cfg_topo{:});
% % ft_plot_lay(lay, cfg_lay{:})
% % title('Spatial filter w')
% %% Source-level ERP
% cfg            = [];
% cfg.keeptrials = 'no';
% avgswitch_LDA      = ft_timelockanalysis(cfg, switchlda);
% avgstay_LDA      = ft_timelockanalysis(cfg, staylda);
% 
% cfg            = [];
% cfg.keeptrials = 'yes';
% tlswitch_LDA       = ft_timelockanalysis(cfg, switchlda);
% tlstay_LDA       = ft_timelockanalysis(cfg, staylda);
% %     cfg                  = [];
% %     cfg.spmversion       = 'spm12';
% %     cfg.method           = 'analytic';
% %     cfg.statistic        = 'ft_statfun_indepsamplesT'; % use the independent samples T-statistic as a measure to
% %     cfg.correctm         = 'fdr';t
% %     cfg.alpha            = 0.05;               % alpha level of the permutation test
% %     cfg.numrandomization = 1000;      % number of draws from the permutation distribution
% %
% %     design                                                                                 = zeros(1,size(tlswitch_LDA.trial,1) + size(tlstay_LDA.trial,1));
% %     design(1,1:size(tlswitch_LDA.trial,1))                                                     = 1;
% %     design(1,(size(tlswitch_LDA.trial,1)+1):(size(tlswitch_LDA.trial,1) + size(tlstay_LDA.trial,1))) = 2;
% %
% %     cfg.design = design;             % design matrix
% %     cfg.ivar   = 1;                   % number or list with indices indicating the independent variable(s)
% %
% %     stat = ft_timelockstatistics(cfg, tlswitch_LDA, tlstay_LDA);
% %
% %     figure;ft_singleplotER([],avgstay_LDA,avgswitch_LDA)
% %     hold on
% %     sig=zeros(length(stat.mask),1);
% %     sig(~stat.mask)=NaN;
% %     plot(stat.time,sig,'* k')
% all_subs_LDAstay{i}=avgstay_LDA;
% all_subs_LDAswitch{i}=avgswitch_LDA;
% end
% 
% 
% % [w,switchlda,C] = LDAbeamformer(P,allswitch,.1);
% % [w,staylda,C] = LDAbeamformer(P,allstay,.1);
% % cfg            = [];
% % cfg.keeptrials = 'no';
% % avgswitch_LDA      = ft_timelockanalysis(cfg, switchlda);
% % avgstay_LDA      = ft_timelockanalysis(cfg, staylda);
% % cfg            = [];
% % cfg.keeptrials = 'yes';
% % tlswitch_LDA       = ft_timelockanalysis(cfg, switchlda);
% % tlstay_LDA       = ft_timelockanalysis(cfg, staylda);
% % all_subs_LDAstay{i}=avgstay_LDA;
% % all_subs_LDAswitch{i}=avgswitch_LDA;
% 
% 

close all
clear all

global ft_default
ft_default.spmversion = 'spm12'; % Force SPM12, SPM8 doesn't go well with mac + 2017b
ft_defaults % This loads the rest of the defaults
files = dir('*.mat');

decode=[];

for i=1:length(files)
load(files(i).name)

% stay_trials   = [events_allBlocks.cuechstay;events_allBlocks.cueenstay];
% switch_trials = [events_allBlocks.cuechswitch;events_allBlocks.cueenswitch];

stay_trials   = [events_allBlocks.cuechstay;events_allBlocks.cueenstay];
switch_trials = [events_allBlocks.cuechswitch;events_allBlocks.cueenswitch];

% ncond_1 = numel(stay_trials);
% ncond_2 = numel(switch_trials);

cfg         = [];
cfg.trials  = stay_trials;
cfg.latency = [-0.25 0.8];
allstay     = ft_selectdata(cfg,all_blocks_clean);

cfg         = [];
cfg.trials  = find(any(~isnan(cell2mat(cellfun(@(isgood)isgood(:),allstay.trial,'uni',0))),1));
allstay     = ft_selectdata(cfg,allstay);

%%%
%trigger art
%%%%

trigger_window = [-0.08 -0.03]; % The window we would like to replace with nans
trigger_window_idx = [nearest(allstay.time{1},trigger_window(1)) nearest(allstay.time{1},trigger_window(2))]; % Find the indices in the time vector corresponding to our window of interest
for j=1:numel(allstay.trial) % Loop through all trials
  allstay.trial{j}(:,trigger_window_idx(1):trigger_window_idx(2))=nan; % Replace the segment of data corresponding to our window of interest with nans
end

% Interpolate nans using cubic interpolation
cfg = [];
cfg.method = 'pchip'; % Here you can specify any method that is supported by interp1: 'nearest','linear','spline','pchip','cubic','v5cubic'
cfg.prewindow = 0.025; % Window prior to segment to use data points for interpolation
cfg.postwindow = 0.025; % Window after segment to use data points for interpolation
allstay = ft_interpolatenan(cfg, allstay); % Clean data

for k=1:length(allstay.trial)
allstay.trial{k} = ft_preproc_smooth(allstay.trial{k}, 7);
end

cfg         = [];
cfg.trials  = switch_trials;
cfg.latency = [-0.25 0.8];
allswitch   = ft_selectdata(cfg,all_blocks_clean);

cfg         = [];
cfg.trials  = find(any(~isnan(cell2mat(cellfun(@(isgood)isgood(:),allswitch.trial,'uni',0))),1));
allswitch     = ft_selectdata(cfg,allswitch);

for j=1:numel(allswitch.trial) % Loop through all trials
  allswitch.trial{j}(:,trigger_window_idx(1):trigger_window_idx(2))=nan; % Replace the segment of data corresponding to our window of interest with nans
end

% Interpolate nans using cubic interpolation
cfg = [];
cfg.method = 'pchip'; % Here you can specify any method that is supported by interp1: 'nearest','linear','spline','pchip','cubic','v5cubic'
cfg.prewindow = 0.025; % Window prior to segment to use data points for interpolation
cfg.postwindow = 0.025; % Window after segment to use data points for interpolation
allswitch = ft_interpolatenan(cfg, allswitch); % Clean data

for k=1:length(allswitch.trial)
allswitch.trial{k} = ft_preproc_smooth(allswitch.trial{k}, 7);
end
% alldata=ft_appenddata([],allstay,allswitch);
% 
% allswitch_ave=ft_timelockanalysis([],allswitch);
% allstay_ave=ft_timelockanalysis([],allstay);
% alldata_ave=ft_timelockanalysis([],alldata);
% 
% cfg=[];
% cfg.operation='subtract';
% cfg.parameter='avg';
% stayvsswitch=ft_math(cfg,allstay_ave,allswitch_ave);
% 
% allswitch_aveall{i} = allswitch_ave;
% allstay_aveall{i} = allstay_ave;
% 
% ival   = [0.15 .25];
% topo_dat = alldata_ave;
% time_idx = (topo_dat.time >= ival(1) & topo_dat.time <= ival(2));
% 
% P        = mean(topo_dat.avg(:, time_idx),2);
% %% LDA beamforming ------------------------------------------------
% % To obtain the LDA beamformer, we need the spatial pattern P and the
% % corresponding continuous or epoched data (NOT the averaged data).
% [w,switchlda,C] = LDAbeamformer(P,allswitch,.1);
% [w,staylda,C] = LDAbeamformer(P,allstay,.1);
% %% Plot covariance matrix, spatial pattern P, and spatial filter w
% % Get channel layout for plotting
% % cfg              = [];
% % cfg.grad         = topo_dat.grad;
% % cfg.skipscale    = 'yes';
% % cfg.skipcomnt    = 'yes';
% % cfg.grad.balance = rmfield(cfg.grad.balance,'reject');
% % lay              = ft_prepare_layout(cfg);
% %
% % % Use the low-level plotting function for plotting arbitrary vectors
% % cfg_topo = {'mask' lay.mask, 'datmask' [], 'interplim' 'mask'};
% % cfg_lay  = { 'point' 'no' 'box' 'no' 'label' 'no' };
% %
% % figure
% % subplot(1,3,1),
% % imagesc(C),title('Covariance matrix')
% % subplot(1,3,2),
% % ft_plot_topo(lay.pos(:,1), lay.pos(:,2), P, cfg_topo{:});
% % ft_plot_lay(lay, cfg_lay{:})
% % title('Spatial pattern P')
% % subplot(1,3,3),
% % ft_plot_topo(lay.pos(:,1), lay.pos(:,2), w, cfg_topo{:});
% % ft_plot_lay(lay, cfg_lay{:})
% % title('Spatial filter w')
% %% Source-level ERP
% cfg            = [];
% cfg.keeptrials = 'no';
% avgswitch_LDA      = ft_timelockanalysis(cfg, switchlda);
% avgstay_LDA      = ft_timelockanalysis(cfg, staylda);
% 
% % cfg            = [];
% % cfg.keeptrials = 'yes';
% % tlswitch_LDA       = ft_timelockanalysis(cfg, switchlda);
% % tlstay_LDA       = ft_timelockanalysis(cfg, staylda);
% %     cfg                  = [];
% %     cfg.spmversion       = 'spm12';
% %     cfg.method           = 'analytic';
% %     cfg.statistic        = 'ft_statfun_indepsamplesT'; % use the independent samples T-statistic as a measure to
% %     cfg.correctm         = 'fdr';t
% %     cfg.alpha            = 0.05;               % alpha level of the permutation test
% %     cfg.numrandomization = 1000;      % number of draws from the permutation distribution
% %
% %     design                                                                                 = zeros(1,size(tlswitch_LDA.trial,1) + size(tlstay_LDA.trial,1));
% %     design(1,1:size(tlswitch_LDA.trial,1))                                                     = 1;
% %     design(1,(size(tlswitch_LDA.trial,1)+1):(size(tlswitch_LDA.trial,1) + size(tlstay_LDA.trial,1))) = 2;
% %
% %     cfg.design = design;             % design matrix
% %     cfg.ivar   = 1;                   % number or list with indices indicating the independent variable(s)
% %
% %     stat = ft_timelockstatistics(cfg, tlswitch_LDA, tlstay_LDA);
% %
% %     figure;ft_singleplotER([],avgstay_LDA,avgswitch_LDA)
% %     hold on
% %     sig=zeros(length(stat.mask),1);
% %     sig(~stat.mask)=NaN;
% %     plot(stat.time,sig,'* k')
% all_subs_LDAstay{i}=avgstay_LDA;
% all_subs_LDAswitch{i}=avgswitch_LDA;


    cfg                 = [] ;
    cfg.method          = 'mvpa';
    cfg.mvpa.classifier = 'lda';
    cfg.mvpa.balance = 'undersample';
    cfg.mvpa.normalise = 'none';
    cfg.mvpa.metric     = {'accuracy', 'auc'};
    cfg.mvpa.k          = 10;
    cfg.mvpa.repeat     = 2;
    cfg.design          = [ones(length(allswitch.trial),1); 2*ones(length(allstay.trial),1)];
    
    stat = ft_timelockstatistics(cfg, allswitch, allstay);
    
    %mv_plot_result(stat.mvpa, allstay.time{1})
    %title(files(i).name)
    
    decode{i} = stat.mvpa;
end


% close all
% all_subs_LDAstay_GM=ft_timelockgrandaverage([],all_subs_LDAstay{:});
% all_subs_LDAswitch_GM=ft_timelockgrandaverage([],all_subs_LDAswitch{:});
% 
% 
% figure;ft_singleplotER([],all_subs_LDAstay_GM,all_subs_LDAswitch_GM)

%figure;plot(allstay.time,mean(decodex))
props = {'-r';'-g';'-m';'-k';'-b';'-y'};

decodex=[];

for i=1:length(decode)
decodex(i,:)=decode{i}.perf{1};
end

[h,p]=ttest(decodex,.5);
p2=FDR(p,0.05);
p(p>p2)=NaN;
p(~isnan(p))=0.55;

figure;shadedErrorBar(allstay.time{1},mean(decodex),...
    std(decodex)/sqrt(size(decode,1)),'lineprops', props{1},'transparent',true,'patchSaturation',0.075)
hold on
plot(allstay.time{1},0.5,'k .')
plot(allstay.time{1},0.5,'. k')
plot(allstay.time{1},p,'* r')
xlim([-0.2 0.8])



% %close all
% figure;plot(allstay.time{1},mean(decodex))
% hold on
% 
% 
% xlim([-0.2 0.8])




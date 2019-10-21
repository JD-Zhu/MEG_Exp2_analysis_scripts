
directories = {'/Volumes/PFS01/ME056/ccsrf_files/Stutterers/CC_files/ANALYSIS_OPERC+TRIANG_MAX_ORI_1mm';...
    '/Volumes/PFS01/ME056/ccsrf_files/Controls_for_Adult_Stuttering_Study/CC_files/ANALYSIS_OPERC+TRIANG_MAX_ORI_5mm'};

for h=1:length(directories)
    cd(directories{h})
    files = dir('ccsrf*mat');
    for i=1:length(files)
        
        D    = spm_eeg_load(files(i).name);
        data = ftraw(D);
        
        cfg        = [];
        cfg.trials = ismember(data.trialinfo,2);
        ignore     = ft_selectdata(cfg,data);
        
        cfg        = [];
        cfg.trials = ismember(data.trialinfo,3);
        goodstop   = ft_selectdata(cfg,data);
        
        cfg        = [];
        cfg.trials = ismember(data.trialinfo,4);
        badstop    = ft_selectdata(cfg,data);
        
        cfg        = [];
        cfg.trials = ismember(data.trialinfo,1);
        go         = ft_selectdata(cfg,data);
        
        cfg     = [];
        allstop = ft_appenddata(cfg,goodstop,badstop,ignore);
        
        cfg         = [];
        av_ignore   = ft_timelockanalysis(cfg,ignore);
        av_goodstop = ft_timelockanalysis(cfg,goodstop);
        av_badstop  = ft_timelockanalysis(cfg,badstop);
        av_allstop  = ft_timelockanalysis(cfg,allstop);
        av_go       = ft_timelockanalysis(cfg,go);
        
        GFP_av_ignore   = ft_globalmeanfield([],av_ignore);
        GFP_av_goodstop = ft_globalmeanfield([],av_goodstop);
        GFP_av_badstop  = ft_globalmeanfield([],av_badstop);
        GFP_av_go       = ft_globalmeanfield([],av_go);
        GFP_av_allstop  = ft_globalmeanfield([],av_allstop);
        
        %Get channel layout for plotting
        cfg           = [];
        cfg.rotate    = 180;
        cfg.grad      = ignore.grad;
        cfg.skipscale = 'yes';
        cfg.skipcomnt = 'yes';
        lay           = ft_prepare_layout(cfg);
        
        cfg           = [];
        cfg.operation = 'subtract';
        cfg.parameter = 'avg';
        %ignorevsgoodstop  = ft_math(cfg,av_ignore,av_goodstop);
        %ignorevsbadstop   = ft_math(cfg,av_ignore,av_badstop);
        %goodstopvsbadstop = ft_math(cfg,av_goodstop,av_badstop);
        inhibition    = ft_math(cfg,av_goodstop,av_ignore); % find the "difference wave" btwn stop & ignore trials
        
        %all_data    = ft_appenddata([],ignore,goodstop,badstop);
        %alldata_avg = ft_timelockanalysis([],all_data);
        
        ival = [0.10 .195;0.19 0.30;0.285 0.39;.445 .7]; % select a few time windows to try, below we will cycle thru these options
        
        % compute the spatial pattern (P1, P2)
        time_idx = (inhibition.time >= ival(1,1) & inhibition.time <= ival(1,2));
        P1       = mean(allstop.avg(:, time_idx),2);
        
        time_idx = (inhibition.time >= ival(2,1) & inhibition.time <= ival(2,2));
        P2       = mean(allstop.avg(:, time_idx),2);
        %% LDA beamforming ------------------------------------------------
        % To obtain the LDA beamformer, we need the spatial pattern P and the
        % corresponding continuous or epoched data (NOT the averaged data).
        [w,t1_goodstoplda,C] = LDAbeamformer(P1,goodstop,.1); % feed in the epoched data for each cond
        [w,t1_badstoplda,C]  = LDAbeamformer(P1,badstop,.1);
        [w,t1_ignorelda,C]   = LDAbeamformer(P1,ignore,.1);
        [w,t1_golda,C]       = LDAbeamformer(P1,go,.1);
        
        [w,t2_goodstoplda,C] = LDAbeamformer(P2,goodstop,.1);
        [w,t2_badstoplda,C]  = LDAbeamformer(P2,badstop,.1);
        [w,t2_ignorelda,C]   = LDAbeamformer(P2,ignore,.1);
        [w,t2_golda,C]       = LDAbeamformer(P2,go,.1);
        
        %% Source-level ERP
        cfg                = [];
        cfg.keeptrials     = 'no';
        t1_avggoodstop_LDA = ft_timelockanalysis(cfg, t1_goodstoplda);
        t1_avgbadstop_LDA  = ft_timelockanalysis(cfg, t1_badstoplda);
        t1_avgignore_LDA   = ft_timelockanalysis(cfg, t1_ignorelda);
        t1_avggo_LDA       = ft_timelockanalysis(cfg, t1_golda);
        
        t2_avggoodstop_LDA = ft_timelockanalysis(cfg, t2_goodstoplda);
        t2_avgbadstop_LDA  = ft_timelockanalysis(cfg, t2_badstoplda);
        t2_avgignore_LDA   = ft_timelockanalysis(cfg, t2_ignorelda);
        t2_avggo_LDA       = ft_timelockanalysis(cfg, t2_golda);
        
        cfg.keeptrials  = 'yes';
        t1_goodstop_LDA = ft_timelockanalysis(cfg, t1_goodstoplda);
        t1_badstop_LDA  = ft_timelockanalysis(cfg, t1_badstoplda);
        t1_ignore_LDA   = ft_timelockanalysis(cfg, t1_ignorelda);
        t1_go_LDA       = ft_timelockanalysis(cfg, t1_golda);
        
        t2_goodstop_LDA = ft_timelockanalysis(cfg, t2_goodstoplda);
        t2_badstop_LDA  = ft_timelockanalysis(cfg, t2_badstoplda);
        t2_ignore_LDA   = ft_timelockanalysis(cfg, t2_ignorelda);
        t2_go_LDA       = ft_timelockanalysis(cfg, t2_golda);
        
        t1_all_goodstop{i} = t1_avggoodstop_LDA;
        t1_all_badstop{i}  = t1_avgbadstop_LDA;
        t1_all_ignore{i}   = t1_avgignore_LDA;
        t1_all_go{i}       = t1_avggo_LDA;
        
        t2_all_goodstop{i} = t2_avggoodstop_LDA;
        t2_all_badstop{i}  = t2_avgbadstop_LDA;
        t2_all_ignore{i}   = t2_avgignore_LDA;
        t2_all_go{i}       = t2_avggo_LDA;
        
        
        t1_all_goodstop_trials{i} = t1_goodstop_LDA;
        t1_all_badstop_trials{i}  = t1_badstop_LDA;
        t1_all_ignore_trials{i}   = t1_ignore_LDA;
        t1_all_go_trials{i}       = t1_go_LDA;
        
        t2_all_goodstop_trials{i} = t2_goodstop_LDA;
        t2_all_badstop_trials{i}  = t2_badstop_LDA;
        t2_all_ignore_trials{i}   = t2_ignore_LDA;
        t2_all_go_trials{i}       = t2_go_LDA;
        
    end
    
    save t1_all_go t1_all_go
    save t2_all_go t2_all_go
    save t1_all_goodstop t1_all_goodstop
    save t2_all_goodstop t2_all_goodstop
    save t2_all_badstop t2_all_badstop
    save t2_all_ignore t2_all_ignore
    save t1_all_ignore t1_all_ignore
    save t1_all_badstop t1_all_badstop
    
    save t1_all_go_trials t1_all_go_trials
    save t2_all_go_trials t2_all_go_trials
    save t1_all_goodstop_trials t1_all_goodstop_trials
    save t2_all_goodstop_trials t2_all_goodstop_trials
    save t1_all_badstop_trials t1_all_badstop_trials
    save t2_all_badstop_trials t2_all_badstop_trials
    save t1_all_ignore_trials t1_all_ignore_trials
    save t2_all_ignore_trials t2_all_ignore_trials
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%snips for group
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plotting

cd(directories{1})
load('t1_all_badstop.mat')
load('t1_all_goodstop.mat')
load('t2_all_goodstop.mat')
load('t1_all_ignore.mat')
load('t2_all_badstop.mat')
load('t2_all_ignore.mat')
load('t1_all_go.mat')
load('t2_all_go.mat')

load('t1_all_badstop_trials.mat')
load('t1_all_goodstop_trials.mat')
load('t2_all_goodstop_trials.mat')
load('t1_all_ignore_trials.mat')
load('t2_all_badstop_trials.mat')
load('t2_all_ignore_trials.mat')
load('t1_all_go_trials.mat')
load('t2_all_go_trials.mat')

Stutt_t1_all_ignore   = t1_all_ignore;
Stutt_t1_all_goodstop = t1_all_goodstop;
Stutt_t2_all_goodstop = t2_all_goodstop;
Stutt_t1_all_badstop  = t1_all_badstop;
Stutt_t2_all_ignore   = t2_all_ignore;
Stutt_t2_all_badstop  = t2_all_badstop;
Stutt_t1_all_go       = t1_all_go;
Stutt_t2_all_go       = t2_all_go;

Stutt_t1_all_ignore_trials   = t1_all_ignore_trials;
Stutt_t1_all_goodstop_trials = t1_all_goodstop_trials;
Stutt_t2_all_goodstop_trials = t2_all_goodstop_trials;
Stutt_t1_all_badstop_trials  = t1_all_badstop_trials;
Stutt_t2_all_ignore_trials   = t2_all_ignore_trials;
Stutt_t2_all_badstop_trials  = t2_all_badstop_trials;
Stutt_t1_all_go_trials       = t1_all_go_trials;
Stutt_t2_all_go_trials       = t2_all_go_trials;

StuttGM_t1_goodstop = ft_timelockgrandaverage([],t1_all_goodstop{:});
StuttGM_t1_badstop  = ft_timelockgrandaverage([],t1_all_badstop{:});
StuttGM_t1_ignore   = ft_timelockgrandaverage([],t1_all_ignore{:});
StuttGM_t1_go       = ft_timelockgrandaverage([],t1_all_go{:});

StuttGM_t2_goodstop = ft_timelockgrandaverage([],t2_all_goodstop{:});
StuttGM_t2_badstop  = ft_timelockgrandaverage([],t2_all_badstop{:});
StuttGM_t2_ignore   = ft_timelockgrandaverage([],t2_all_ignore{:});
StuttGM_t2_go       = ft_timelockgrandaverage([],t2_all_go{:});

cd(directories{2})
load('t1_all_badstop.mat')
load('t2_all_badstop.mat')
load('t1_all_goodstop.mat')
load('t2_all_goodstop.mat')
load('t1_all_ignore.mat')
load('t2_all_ignore.mat')
load('t1_all_go.mat')
load('t2_all_go.mat')

load('t1_all_badstop_trials.mat')
load('t1_all_goodstop_trials.mat')
load('t2_all_goodstop_trials.mat')
load('t1_all_ignore_trials.mat')
load('t2_all_badstop_trials.mat')
load('t2_all_ignore_trials.mat')
load('t1_all_go_trials.mat')
load('t2_all_go_trials.mat')

Ctrl_t1_all_ignore   = t1_all_ignore;
Ctrl_t1_all_goodstop = t1_all_goodstop;
Ctrl_t2_all_goodstop = t2_all_goodstop;
Ctrl_t1_all_badstop  = t1_all_badstop;
Ctrl_t2_all_ignore   = t2_all_ignore;
Ctrl_t2_all_badstop  = t2_all_badstop;
Ctrl_t1_all_go       = t1_all_go;
Ctrl_t2_all_go       = t2_all_go;

Ctrl_t1_all_ignore_trials   = t1_all_ignore_trials;
Ctrl_t1_all_goodstop_trials = t1_all_goodstop_trials;
Ctrl_t2_all_goodstop_trials = t2_all_goodstop_trials;
Ctrl_t1_all_badstop_trials  = t1_all_badstop_trials;
Ctrl_t2_all_ignore_trials   = t2_all_ignore_trials;
Ctrl_t2_all_badstop_trials  = t2_all_badstop_trials;
Ctrl_t1_all_go_trials       = t1_all_go_trials;
Ctrl_t2_all_go_trials       = t2_all_go_trials;

CtrlGM_t1_goodstop = ft_timelockgrandaverage([],t1_all_goodstop{:});
CtrlGM_t1_badstop  = ft_timelockgrandaverage([],t1_all_badstop{:});
CtrlGM_t1_ignore   = ft_timelockgrandaverage([],t1_all_ignore{:});
CtrlGM_t1_go       = ft_timelockgrandaverage([],t1_all_go{:});

CtrlGM_t2_goodstop = ft_timelockgrandaverage([],t2_all_goodstop{:});
CtrlGM_t2_badstop  = ft_timelockgrandaverage([],t2_all_badstop{:});
CtrlGM_t2_ignore   = ft_timelockgrandaverage([],t2_all_ignore{:});
CtrlGM_t2_go       = ft_timelockgrandaverage([],t2_all_go{:});

figure;
subplot(2,4,1)
ft_singleplotER([],StuttGM_t1_goodstop,CtrlGM_t1_goodstop)
ylim([-0.5 1.5])
subplot(2,4,2)
ft_singleplotER([],StuttGM_t1_badstop,CtrlGM_t1_badstop)
ylim([-0.5 1.5])
subplot(2,4,3)
ft_singleplotER([],StuttGM_t1_ignore,CtrlGM_t1_ignore)
ylim([-0.5 1.5])
subplot(2,4,4)
ft_singleplotER([],StuttGM_t1_go,CtrlGM_t1_go)
ylim([-0.5 1.5])

subplot(2,4,5)
ft_singleplotER([],StuttGM_t2_goodstop,CtrlGM_t2_goodstop)
ylim([-0.5 1.5])
subplot(2,4,6)
ft_singleplotER([],StuttGM_t2_badstop,CtrlGM_t2_badstop)
ylim([-0.5 1.5])
subplot(2,4,7)
ft_singleplotER([],StuttGM_t2_ignore,CtrlGM_t2_ignore)
ylim([-0.5 1.5])
subplot(2,4,8)
ft_singleplotER([],StuttGM_t2_go,CtrlGM_t2_go)
ylim([-0.5 1.5])

legend('stutt','ctrl')

cfg                  = [];
cfg.spmversion       = 'spm12';
cfg.method           = 'analytic';
cfg.statistic        = 'ft_statfun_indepsamplesT'; % use the independent samples T-statistic as a measure to
cfg.correctm         = 'fdr';
cfg.alpha            = 0.05;               % alpha level of the permutation test
cfg.numrandomization = 1000;      % number of draws from the permutation distribution

design                                                                  = zeros(1,size(t1_all_go,2) + size(t1_all_go,2));
design(1,1:size(t1_all_go,2))                                           = 1;
design(1,(size(t1_all_go,2)+1):(size(t1_all_go,2) + size(t1_all_go,2))) = 2;

cfg.design = design;             % design matrix
cfg.ivar   = 1;                   % number or list with indices indicating the independent variable(s)

stat = ft_timelockstatistics(cfg, Stutt_t2_all_badstop{:}, Ctrl_t2_all_badstop{:});

figure;ft_singleplotER([],StuttGM_t2_badstop,CtrlGM_t2_badstop)
hold on
sig                      = zeros(length(stat.prob),1);
sig(fdr(stat.prob)<0.05) = 1;
sig(~sig)                = NaN;
plot(Ctrl_t1_all_go{1}.time,sig,'* k')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Ctrl_t1_all_go_tf  = [];
% Stutt_t1_all_go_tf = [];
% %
% for i =1:length(Ctrl_t1_all_go_trials)
%     
%     cfg                  = [];
%     cfg.output           = 'pow';
%     cfg.channel          = 'all';
%     cfg.method           = 'mtmconvol';
%     cfg.taper            = 'hanning';
%     cfg.foi              = 2:0.5:40;                         % analysis 5 to 30 Hz in steps of 0.5 Hz
%     cfg.t_ftimwin        = 5./cfg.foi;
%     cfg.tapsmofrq        = 0.4 * cfg.foi;
%     cfg.toi              = -0.5:0.05:1.0;                  % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
%     cfg.pad              = 'nextpow2';
%     go_ft_ctrl           = ft_freqanalysis(cfg, Ctrl_t1_all_go_trials{i});
%     Ctrl_t1_all_go_tf{i} = go_ft_ctrl;
% end
% 
% for i =1:length(Stutt_t1_all_go_trials)
%     
%     cfg                   = [];
%     cfg.output            = 'pow';
%     cfg.channel           = 'all';
%     cfg.method            = 'mtmconvol';
%     cfg.taper             = 'hanning';
%     cfg.foi               = 2:0.5:40;                         % analysis 5 to 30 Hz in steps of 0.5 Hz
%     cfg.t_ftimwin         = 5./cfg.foi;
%     cfg.tapsmofrq         = 0.4 * cfg.foi;
%     cfg.toi               = -0.5:0.05:1.0;                  % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
%     cfg.pad               = 'nextpow2';
%     go_ft_stutt           = ft_freqanalysis(cfg, Stutt_t1_all_go_trials{i});
%     Stutt_t1_all_go_tf{i} = go_ft_stutt;
% end
% 
% 
%   cfg = [];
%   cfg.spmversion       = 'spm12';
%   cfg.channel          =  'LDAsource';
%   cfg.frequency        = [5 30];
%   cfg.neighbourdist    = 4;
%   cfg.latency          = [0 .5];
%   cfg.avgovertime      = 'no';
%   cfg.avgoverfreq      ='no';
%   cfg.avgoverchan      = 'no';
%  
%   cfg.clusteralpha     = 0.05;
%   cfg.statistic        = 'ft_statfun_indepsamplesT'; % use the independent samples T-statistic as a measure to
%   cfg.numrandomization = 500;
%   cfg.correctm         =   'none';
%   cfg.method           = 'analytic';
% 
%   
% design                                                                  = zeros(1,size(Ctrl_t1_all_go_tf,2) + size(Ctrl_t1_all_go_tf,2));
% design(1,1:size(Ctrl_t1_all_go_tf,2))                                           = 1;
% design(1,(size(Ctrl_t1_all_go_tf,2)+1):(size(Ctrl_t1_all_go_tf,2) + size(Ctrl_t1_all_go_tf,2))) = 2;
% 
% cfg.design = design;             % design matrix
% cfg.ivar   = 1;                   % number or list with indices indicating the independent variable(s)
% 
% stat = ft_freqstatistics(cfg, Ctrl_t1_all_go_tf{:}, Stutt_t1_all_go_tf{:});
% %
% Ctrl_t1_all_go_tf_bc  = [];
% Stutt_t1_all_go_tf_bc = [];
% cfg                   = [];
% cfg.baseline          = [-0.2 0];
% cfg.baselinetype      = 'relchange';
% 
% for i=1:length(Ctrl_t1_all_go_tf)
%     Ctrl_t1_all_go_tf_bc{i} = ft_freqbaseline(cfg,Ctrl_t1_all_go_tf{i});
% end
% 
% for i=1:length(Stutt_t1_all_go_tf)
%     Stutt_t1_all_go_tf_bc{i} = ft_freqbaseline(cfg,Stutt_t1_all_go_tf{i});
% end
% 
% Ctrl_t1_all_go_tf_GM  = ft_freqgrandaverage([],Ctrl_t1_all_go_tf_bc{:});
% Stutt_t1_all_go_tf_GM = ft_freqgrandaverage([],Stutt_t1_all_go_tf_bc{:});
% 
% figure;
% cfg      = [];
% cfg.xlim = [-0.2 0.75];
% cfg.ylim = [5 30];
% cfg.zlim = [-.25 .25];
% subplot(2,1,1)
% ft_singleplotTFR(cfg,Ctrl_t1_all_go_tf_GM)
% subplot(2,1,2)
% ft_singleplotTFR(cfg,Stutt_t1_all_go_tf_GM)

%
%
% go_ft_ctrl         = ft_freqanalysis(cfg, Ctrl_t1_all_goodstop{1});
%
%
% badstop          = ft_freqanalysis(cfg, Ctrl_t1_all_badstop{1});
% % go          = ft_freqanalysis(cfg, Ctrl_t1_all_go{1});
% % %
% % cfg              = [];
% cfg.baseline     = [-2.0 -1.1];
% cfg.baselinetype = 'relchange';
% cfg.zlim         = [-0.3 0];
% cfg.xlim         = [-0.5 1];
% cfg.ylim         = [2 30];
% cfg.channel      = 'all'; % top figure
%
% figure
% subplot(2,1,1)
% ft_singleplotTFR(cfg, goodstop);
% title('Goodstop')
% subplot(2,1,2)
% ft_singleplotTFR(cfg, badstop);
% title('Badstop')

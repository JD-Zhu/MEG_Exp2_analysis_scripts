%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% stats_ERF.m
%
% Author: Judy Zhu (github.com/JD-Zhu)
%
% Grand average & statistical analysis on sensor-space ERFs.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
%clear all;


% = Settings =
% Please adjust as required:

% perform channel repair on each subject's ERF?
CHANNEL_REPAIR = false; % only need to do the repair once, we'll save repaired erf in the folder below
repaired_erf_folder = 'channelrepaired\\'; % need to create this folder first

% cfg.avgovertime setting in cluster-based permutation test
AVGOVERTIME = true;
TIME_WINDOW_TO_AVG = [0.060 0.110]; % must set this var, if AVGOVERTIME is set to true


% SELECT which set of single-subject ERFs to use
run_name = 'TSPCA10000_3'; % this should be a folder name inside the "Results_ERF" folder



%%
% run the #define section
global ResultsFolder; % all subjects' erf data are stored here
global filename_suffix; % erf results file suffix

global eventnames_real; global colours_and_lineTypes; 
global colours; global lineTypes;
global PLOT_XLIM; global ERF_BASELINE;
common();


% initialise allSubjects_erf (each field holds all subjects' erf in that condition)
allSubjects_erf.NatStay = {};
allSubjects_erf.NatSwitch = {};
allSubjects_erf.NatSingle = {};
allSubjects_erf.ArtStay = {};
allSubjects_erf.ArtSwitch = {};
allSubjects_erf.ArtSingle = {};
allSubjects_erf.BiStay = {};
allSubjects_erf.BiSwitch = {};
allSubjects_erf.BiSingle = {};


%% Read data

% find all .mat files in ResultsFolder_thisrun
ResultsFolder_thisrun = [ResultsFolder run_name '\\']; % ERF results for all subjects
files = dir([ResultsFolder_thisrun '*_erf' filename_suffix '.mat']);

% each cycle reads in one '.mat' file (ie. one subject's erf data)
for i = 1:length(files)
    filename = [ResultsFolder_thisrun files(i).name];
    load(filename);
        
    for j = 1:length(eventnames_real) % 9 conds (if collapsed across langs) or 18 conds (if not collapsed)
        % perform channel repair if needed
        if (CHANNEL_REPAIR == true)
            load('neighbours.mat');
            load('all_labels.mat');
            erf_clean.(eventnames_real{j}) = repair_bad_channels(erf_clean.(eventnames_real{j}), neighbours, all_labels);
        end
        % add to allsubjects matrix
        allSubjects_erf.(eventnames_real{j}) = [allSubjects_erf.(eventnames_real{j}) erf_clean.(eventnames_real{j})];
    end
    
    % save the new erf after channel repair
    if (CHANNEL_REPAIR == true)
        save([ResultsFolder_thisrun repaired_erf_folder files(i).name], 'SubjectFolder', 'erf_clean');
    end
end


%% Descriptives
% http://www.fieldtriptoolbox.org/tutorial/cluster_permutation_timelock#within-subjects_experiments

fprintf('\n= COMPUTING & PLOTTING CROSS-SUBJECT AVERAGES =\n');

% CALCULATE the grand average (across all subjects) for each condition
cfg = [];
cfg.channel   = {'all', '-AG101', '-AG122', '-AG007', '-AG103'}; % remove noisy sensors (MEG Exp2):
                                                                % ch100 (AG101) is always noisy -> Remove for all ptps!
                                                                % ch006 & ch102 also shows the same noise occasionally.
                                                                % ch121 (AG122) tends to show square noise.
cfg.latency   = 'all';
cfg.parameter = 'avg';
for j = 1:length(eventnames_real)
    cfg.keepindividual = 'no'; % average across subjects
    GA_erf.(eventnames_real{j}) = ft_timelockgrandaverage(cfg, allSubjects_erf.(eventnames_real{j}){:});  

    cfg.keepindividual = 'yes'; % do not average across subjects, keep the data for each individual subject
    GA_indi.(eventnames_real{j}) = ft_timelockgrandaverage(cfg, allSubjects_erf.(eventnames_real{j}){:}); 

    % "{:}" means to use data from all elements of the variable
end

% save the GA files
GA_output_file = [ResultsFolder_thisrun 'GA_avg.mat'];
if (exist(GA_output_file, 'file') ~= 2) 
    save(GA_output_file, 'GA_erf');
end
GA_output_file = [ResultsFolder_thisrun 'GA_individuals.mat'];
if (exist(GA_output_file, 'file') ~= 2) 
    save(GA_output_file, 'GA_indi');
end

% multiplot
load('lay.mat');
        
cfg              = [];
cfg.showlabels   = 'yes';
cfg.fontsize     = 6;
cfg.layout       = lay;
cfg.baseline     = ERF_BASELINE;
cfg.baselinetype = 'absolute';
cfg.graphcolor   = cell2mat(colours); 
cfg.linestyle    = lineTypes;
cfg.xlim         = PLOT_XLIM;

figure('Name','ft_multiplotER: GA_erf (9 conds)');
% convert struct to cell array, then you can feed it in as 'varargin'
cellarray = struct2cell(GA_erf);
ft_multiplotER(cfg, cellarray{:});
% specify the legends manually (otherwise it will display incorrectly)
lines = findall(gcf, 'Type','line');
lines = lines([29 26 23 20 17 14 11 8 5]); % grab the correct lines (this is complicated because many lines are plotted in ft_multiplot)
legend(lines, eventnames_real);


% CALCULATE global averages across all sensors (i.e. GFP = global field power)
cfg        = [];
cfg.method = 'power';
%cfg.channel = {'AG017', 'AG018', 'AG019', 'AG022', 'AG023', 'AG025', 'AG029', 'AG063', 'AG064', 'AG143'}; % 10 sig channels in cluster
cfg.channel   = {'all', '-AG101', '-AG122', '-AG007', '-AG103'}; % remove noisy sensors (see above)
for j = 1:length(eventnames_real)
    GA_erf_GFP.(eventnames_real{j}) = ft_globalmeanfield(cfg, GA_erf.(eventnames_real{j}));
end

% plot GFP
figure('Name','GFP_all_subjects'); hold on
for j = 1:length(eventnames_real)
    if colours_and_lineTypes % use a combination of colours and line types to distinguish conds
        plot(GA_erf_GFP.(eventnames_real{j}).time, GA_erf_GFP.(eventnames_real{j}).avg, 'color', colours{j}, 'LineStyle', lineTypes{j});
    else % just use colours
        plot(GA_erf_GFP.(eventnames_real{j}).time, GA_erf_GFP.(eventnames_real{j}).avg, 'color', colours(j,:));
    end
    xlim(PLOT_XLIM);
end
legend(eventnames_real);


% average across all 4 conds (for selecting windows for peaks)
averageAcrossConds = GA_erf_GFP.NatStay;
averageAcrossConds.avg = (GA_erf_GFP.NatStay.avg + GA_erf_GFP.NatSwitch.avg + GA_erf_GFP.NatSingle.avg ...
                        + GA_erf_GFP.ArtStay.avg + GA_erf_GFP.ArtSwitch.avg + GA_erf_GFP.ArtSingle.avg ...
                        + GA_erf_GFP.BiStay.avg + GA_erf_GFP.BiSwitch.avg + GA_erf_GFP.BiSingle.avg) / 9;

figure('Name','GFP_all_subjects - Averaged across all conds'); 
plot(averageAcrossConds.time, averageAcrossConds.avg); 
xlim([-0.2 0.8]);


%% Statistical analysis

data = allSubjects_erf; % make an easy name

% Optional: select which subjects to (not) use
for j = 1:length(eventnames_real)
    %data.(eventnames_real{j})([1 11 16 19 20 21 14]) = [];
end


fprintf('\n= STATS: CLUSTER-BASED PERMUTATION TESTS =\n');

cfg = [];
cfg.channel   = {'all', '-AG101', '-AG122', '-AG007', '-AG103'}; % remove noisy sensors (see above)
load('neighbours.mat'); % this is the sensor layout - it's the same for all subjects (even same across experiments). So just prepare once & save, then load here
cfg.neighbours = neighbours;  % same as defined for the between-trials experiment

% can choose diff time windows to analyse for cue epochs & target epochs
% (these will be fed into cfg.latency accordingly)
if (AVGOVERTIME)
    latency_cue = TIME_WINDOW_TO_AVG; % time range to average over
    cfg.avgovertime = 'yes'; % if yes, this will average over the entire time window chosen in cfg.latency 
                            % (useful when you want to look at a particular component, e.g. to look at M100,
                            % cfg.latency = [0.08 0.12]; cfg.avgovertime = 'yes'; )
else % autoly detect temporal cluster
    latency_cue = [-0.1 0.75]; % time interval over which the experimental 
                               % conditions must be compared (in seconds)
    cfg.avgovertime = 'no';
end

cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT'; %cfg.statistic = 'ft_statfun_indepsamplesT'; OR 'ft_statfun_depsamplesFmultivariate';
cfg.correctm = 'cluster'; %'no'; % its common in MEG studies to run uncorrected at cfg.alpha = 0.001
cfg.clusteralpha = 0.05; % threshold for selecting candidate samples to form clusters
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 2; % minimum number of neighbourhood channels required to be significant 
                   % in order to form a cluster 
                   % (default: 0, ie. each single channel can be considered a cluster).
                   % 4 or 5 is a good choice; 2 is too few coz it's even below
                   % the resolution of the sensor layout (i.e. 2 adjacent sensors might
                   % really be measuring the same thing, so ofc they are both sig)

cfg.tail = 0;
cfg.clustertail = 0; % 2 tailed test
cfg.alpha = 0.1; %0.001  % threshold for cluster-level statistics (any cluster with a p-value lower than this will be reported as sig - an entry of '1' in .mask field)
cfg.correcttail = 'prob'; % correct for 2-tailedness
cfg.numrandomization = 500; % Rule of thumb: use 500, and double this number if it turns out 
    % that the p-value differs from the chosen alpha (e.g. 0.05) by less than 0.02

numSubjects = length(data.(eventnames_real{1})); % check how many subjects we are including
within_design_2x2 = zeros(2,2*numSubjects);
within_design_2x2(1,:) = repmat(1:numSubjects,1,2);
within_design_2x2(2,1:numSubjects) = 1;
within_design_2x2(2,numSubjects+1:2*numSubjects) = 2;

cfg.design = within_design_2x2;
cfg.uvar  = 1; % row of design matrix that contains unit variable (in this case: subjects)
cfg.ivar  = 2; % row of design matrix that contains independent variable (i.e. the conditions)

cfg.latency = latency_cue;


% Run the statistical tests
[Bi_sw] = ft_timelockstatistics(cfg, data.BiStay{:}, data.BiSwitch{:}); 
length(find(Bi_sw.mask))
%%
% Switch cost in each context
[Nat_sw] = ft_timelockstatistics(cfg, data.NatStay{:}, data.NatSwitch{:}); %allSubj_cue_ch_switchCost{:}, allSubj_cue_en_switchCost{:});
[Art_sw] = ft_timelockstatistics(cfg, data.ArtStay{:}, data.ArtSwitch{:});
[Bi_sw] = ft_timelockstatistics(cfg, data.BiStay{:}, data.BiSwitch{:}); 

% Mixing cost in each context
[Nat_mix] = ft_timelockstatistics(cfg, data.NatSingle{:}, data.NatStay{:});
[Art_mix] = ft_timelockstatistics(cfg, data.ArtSingle{:}, data.ArtStay{:});
[Bi_mix] = ft_timelockstatistics(cfg, data.BiSingle{:}, data.BiStay{:}); 

length(find(Nat_sw.mask))
length(find(Art_sw.mask))
length(find(Bi_sw.mask))
length(find(Nat_mix.mask))
length(find(Art_mix.mask))
length(find(Bi_mix.mask))

%save([ResultsFolder_thisrun 'stats_minnbchan' mat2str(cfg.minnbchan) '.mat'], 'Nat_sw', 'Art_sw', 'Bi_sw', 'Nat_mix', 'Art_mix', 'Bi_mix');

% interaction in switch cost
cfg.minnbchan = 2;
[timelock1, timelock2] = combine_conds_for_T_test('fieldtrip', 'interaction', data.NatStay, data.NatSwitch, data.BiStay, data.BiSwitch);
[interaction_sw] = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:}); %allSubj_cue_ch_switchCost{:}, allSubj_cue_en_switchCost{:});

length(find(interaction_sw.mask))

% interaction in mixing cost
cfg.minnbchan = 2;
[timelock1, timelock2] = combine_conds_for_T_test('fieldtrip', 'interaction', data.NatStay, data.NatSingle, data.BiStay, data.BiSingle);
[interaction_mix] = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:}); %allSubj_cue_ch_switchCost{:}, allSubj_cue_en_switchCost{:});

length(find(interaction_mix.mask))


%% Below are from MEG Exp 1
% Interaction (i.e. calc sw$ in each lang, then test the 2 sw$)
% http://www.fieldtriptoolbox.org/faq/how_can_i_test_an_interaction_effect_using_cluster-based_permutation_tests
%{
for i = 1:numSubjects
    allSubj_cue_ch_switchCost{i} = allSubjects_erf.cuechstay{i}; % difference btwn chstay & chswitch (calc'd on next line)
    allSubj_cue_ch_switchCost{i}.avg = allSubjects_erf.cuechswitch{i}.avg - allSubjects_erf.cuechstay{i}.avg; % chswitch - chstay
    allSubj_cue_en_switchCost{i} = allSubjects_erf.cueenstay{i};
    allSubj_cue_en_switchCost{i}.avg = allSubjects_erf.cueenswitch{i}.avg - allSubjects_erf.cueenstay{i}.avg; % enswitch - enstay
    allSubj_target_ch_switchCost{i} = allSubjects_erf.targetchstay{i};
    allSubj_target_ch_switchCost{i}.avg = allSubjects_erf.targetchswitch{i}.avg - allSubjects_erf.targetchstay{i}.avg; % chswitch - chstay
    allSubj_target_en_switchCost{i} = allSubjects_erf.targetenstay{i};
    allSubj_target_en_switchCost{i}.avg = allSubjects_erf.targetenswitch{i}.avg - allSubjects_erf.targetenstay{i}.avg; % enswitch - enstay
end
%}
% manual calculation above is now replaced by combine_conds_for_T_Test()
fprintf('\nCUE window -> Testing lang x ttype interaction:\n');
[timelock1, timelock2] = combine_conds_for_T_test('fieldtrip', 'interaction', data.cuechstay, data.cuechswitch, data.cueenstay, data.cueenswitch);
cfg.latency = latency_cue;  
[cue_interaction] = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:}); %allSubj_cue_ch_switchCost{:}, allSubj_cue_en_switchCost{:});
fprintf('\nTARGET window -> Testing lang x ttype interaction:\n');
[timelock1, timelock2] = combine_conds_for_T_test('fieldtrip', 'interaction', data.targetchstay, data.targetchswitch, data.targetenstay, data.targetenswitch); %'2-1 vs 4-3');
cfg.latency = latency_target;  
[target_interaction] = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:}); %allSubj_target_ch_switchCost{:}, allSubj_target_en_switchCost{:});

% Main effect of lang (collapse across stay-switch)
%{
for i = 1:numSubjects
    allSubj_cue_ch{i} = allSubjects_erf.cuechstay{i};
    allSubj_cue_ch{i}.avg = (allSubjects_erf.cuechstay{i}.avg + allSubjects_erf.cuechswitch{i}.avg) / 2; % cue_ch_all.avg = (cuechstay.avg + cuechsw.avg) / 2
    allSubj_cue_en{i} = allSubjects_erf.cueenstay{i};
    allSubj_cue_en{i}.avg = (allSubjects_erf.cueenstay{i}.avg + allSubjects_erf.cueenswitch{i}.avg) / 2;
    allSubj_target_ch{i} = allSubjects_erf.targetchstay{i};
    allSubj_target_ch{i}.avg = (allSubjects_erf.targetchstay{i}.avg + allSubjects_erf.targetchswitch{i}.avg) / 2; % target_ch_all.avg = (targetchstay.avg + targetchsw.avg) / 2
    allSubj_target_en{i} = allSubjects_erf.targetenstay{i};
    allSubj_target_en{i}.avg = (allSubjects_erf.targetenstay{i}.avg + allSubjects_erf.targetenswitch{i}.avg) / 2;    
end
%}
fprintf('\nCUE window -> Main effect of lang:\n');
[timelock1, timelock2] = combine_conds_for_T_test('fieldtrip', 'main_12vs34', data.cuechstay, data.cuechswitch, data.cueenstay, data.cueenswitch);
cfg.latency = latency_cue; 
[cue_lang] = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:}); %allSubj_cue_ch{:}, allSubj_cue_en{:});
fprintf('\nTARGET window -> Main effect of lang:\n');
[timelock1, timelock2] = combine_conds_for_T_test('fieldtrip', 'main_12vs34', data.targetchstay, data.targetchswitch, data.targetenstay, data.targetenswitch); %'2-1 vs 4-3');
cfg.latency = latency_target;
[target_lang] = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:}); %allSubj_target_ch{:}, allSubj_target_en{:});

% Main effect of switch (collapse across langs)
%{
for i = 1:numSubjects
    allSubj_cue_stay{i} = allSubjects_erf.cuechstay{i};
    allSubj_cue_stay{i}.avg = (allSubjects_erf.cuechstay{i}.avg + allSubjects_erf.cueenstay{i}.avg) / 2;
    allSubj_cue_switch{i} = allSubjects_erf.cuechswitch{i};
    allSubj_cue_switch{i}.avg = (allSubjects_erf.cuechswitch{i}.avg + allSubjects_erf.cueenswitch{i}.avg) / 2;
    allSubj_target_stay{i} = allSubjects_erf.targetchstay{i};
    allSubj_target_stay{i}.avg = (allSubjects_erf.targetchstay{i}.avg + allSubjects_erf.targetenstay{i}.avg) / 2;
    allSubj_target_switch{i} = allSubjects_erf.targetchswitch{i};
    allSubj_target_switch{i}.avg = (allSubjects_erf.targetchswitch{i}.avg + allSubjects_erf.targetenswitch{i}.avg) / 2;
end
%}
fprintf('\nCUE window -> Main effect of ttype:\n');
[timelock1, timelock2] = combine_conds_for_T_test('fieldtrip', 'main_13vs24', data.cuechstay, data.cuechswitch, data.cueenstay, data.cueenswitch);
cfg.latency = latency_cue; 
[cue_ttype] = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:}); %allSubj_cue_stay{:}, allSubj_cue_switch{:});
fprintf('\nTARGET window -> Main effect of ttype:\n');
[timelock1, timelock2] = combine_conds_for_T_test('fieldtrip', 'main_13vs24', data.targetchstay, data.targetchswitch, data.targetenstay, data.targetenswitch); %'2-1 vs 4-3');
cfg.latency = latency_target; 
[target_ttype] = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:}); %allSubj_target_stay{:}, allSubj_target_switch{:});

% check for effects by searching the .mask field
effect_cue_interaction = length(find(cue_interaction.mask)) % if not 0, then we have an effect here
effect_cue_lang = length(find(cue_lang.mask))
effect_cue_ttype = length(find(cue_ttype.mask))
effect_target_interaction = length(find(target_interaction.mask))
effect_target_lang = length(find(target_lang.mask))
effect_target_ttype = length(find(target_ttype.mask))

%save([ResultsFolder_thisrun 'stats.mat'], 'cue_interaction', 'cue_lang', 'cue_ttype', 'target_interaction', 'target_lang', 'target_ttype');


%% Plotting: use ft_clusterplot & ft_topoplot

load([ResultsFolder_thisrun 'stats.mat']);
load('lay.mat');
load([ResultsFolder_thisrun 'GA_erf_allConditions.mat']); % only required if using ft_topoplot

% use a nice-looking colourmap
ft_hastoolbox('brewermap', 1); % ensure this toolbox is on the path
cmap = colormap(flipud(brewermap(64, 'RdBu')));

% select which comparison to plot
stat = target_lang; % here we plot the only effect that seems to survive correction (at minnbchan = 0)
                  % to explore where (both in terms of time & location) the effect might have possibly
                  % occurred
                  % [TODO] then we can define more precise time window &
                  % set avgovertime = 'yes', which should give us more
                  % sensitivity, and allow us to increase the minnbchan to
                  % a reasonable number: 2 (ft tutorial) or 4 (Paul)

%% ft_clusterplot (based on t-values)

% too much warning, can't see the console output
ft_warning off 'FieldTrip:ft_clusterplot:ft_topoplotTFR:topoplot_common:ft_selectdata:getdimord:warning_dimord_could_not_be_determined:line681'

cfg = [];
%cfg.zlim = [-5 5]; % set scaling (range of t-values) (usually using automatic is ok) 
cfg.highlightcolorpos = [1 1 1]; % white for pos clusters
cfg.highlightcolorneg = [255/255 192/255 203/255]; % pink for neg clusters
cfg.alpha = 0.1; % any clusters with a p-value below this threshold will be plotted
cfg.layout = lay;
cfg.colormap = cmap;

% turn on the following lines if you are after one particular subplot
%cfg.subplotsize = [1 1];
%cfg.colorbar = 'yes'; % shows the scaling

ft_clusterplot(cfg, stat);


%% ft_topoplot (based on actual erf amplitude) 
%{
% first, define the 2 conds to be compared (this time using cross-subject averages, i.e. GA)
% here we look at main effect of ttype in cue window, so we collapse across langs
GA_cue_stay = GA_erf.cuechstay;
GA_cue_stay.avg = (GA_erf.cuechstay.avg + GA_erf.cueenstay.avg) / 2;
GA_cue_switch = GA_erf.cuechswitch;
GA_cue_switch.avg = (GA_erf.cuechswitch.avg + GA_erf.cueenswitch.avg) / 2;

% then, calc the diff btwn the 2 conds
cfg  = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
GA_cue_stayvsswitch = ft_math(cfg, GA_cue_stay, GA_cue_switch);


% define parameters for plotting
start_time = stat.cfg.latency(1); % get the time window specified earlier in stat analysis
end_time = stat.cfg.latency(end);
timestep = 0.05; %(end_time - start_time) / 15; % length of time interval you want in each subplot (in seconds); 
                                                % alt: specify how many subplots you want (e.g. 15)
sampling_rate = 200; % we downsampled to 200Hz
sample_count = length(stat.time); % number of samples in MEG data (in the ERF time window)
j = [start_time : timestep : end_time];   % define the time interval (in seconds) for each subplot
m = [1 : timestep*sampling_rate : sample_count];  % corresponding sample indices in MEG data

% ensure stat.cfg.alpha (the alpha level we specified earlier in ft_timelockstatistics) still exists
if ~isfield(stat.cfg,'alpha'); stat.cfg.alpha = 0.05; end; % if not, set it (new version corrects for 2-tailedness, so no need to use 0.025)

%{
if (length(stat.posclusters) == 0) % if no clusters were found at all, code below will throw error
    % so create a fake one (just to allow code below to run w/o error)
    stat.posclusters(1).prob = 1; 
    stat.posclusters(1).clusterstat = -9; 
    stat.posclusters(1).stddev = 0; 
    stat.posclusters(1).cirange = 0;
end
if (length(stat.negclusters) == 0) % do the same for neg clusters
    stat.negclusters(1).prob = 1; 
    stat.negclusters(1).clusterstat = -9; 
    stat.negclusters(1).stddev = 0; 
    stat.negclusters(1).cirange = 0;
end
%}

% get all p-values associated with the clusters
pos_cluster_pvals = [stat.posclusters(:).prob];
neg_cluster_pvals = [stat.negclusters(:).prob];
% find which clusters are significant, outputting their indices as held in stat.posclusters
pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha);
neg_signif_clust = find(neg_cluster_pvals < stat.cfg.alpha);
% make a boolean matrix of which (channel,time)-pairs are part of a significant cluster
pos = ismember(stat.posclusterslabelmat, pos_signif_clust); % I think stat.mask is simply combining pos & neg 
neg = ismember(stat.negclusterslabelmat, neg_signif_clust); % (i.e. stat.mask == pos | neg)

% Ensure the channels have the same order in the grand average and in the statistical output
% This might not be the case, because ft_math might shuffle the order  
[i1,i2] = match_str(GA_cue_stayvsswitch.label, stat.label);
% i1 holds a list of channel numbers in the grand averages
% i2 holds a list of channel numbers in the stat output

figure;  
for k = 1:length(j)-1; % create one subplot for each time interval
     subplot(3,5,k); % 3 * 5 = 15 subplots 
     
     cfg = [];   
     cfg.xlim=[j(k) j(k+1)];   
     %cfg.zlim = [-5e-14 5e-14];  % set scaling (usually using automatic is ok) 
     pos_int = zeros(numel(GA_cue_stayvsswitch.label),1); % initialise the arrays with 0s
     neg_int = zeros(numel(GA_cue_stayvsswitch.label),1);
     pos_int(i1) = all(pos(i2, m(k):m(k+1)), 2); % if a channel maintains significance thruout this time interval, then
     neg_int(i1) = all(neg(i2, m(k):m(k+1)), 2); % we set this channel to 1 (to be highlighted)
     % not sure why it has to "maintain significance"; here I try with only requiring sig for half of time pts in this interval
     a = neg(i2, m(k):m(k+1));
     neg_int(i1) = sum(a, 2) > size(a, 2) / 2;
     
     sig_channels = find(pos_int | neg_int); % get indices of all significant channels
     if length(sig_channels) ~= 0 % if any sig channels found, report which channels these are
         fprintf(['In time interval [' num2str(cfg.xlim) '], these channels were significant:\n']);
         stat.label(sig_channels)
     end
     cfg.highlight = 'on';
     cfg.highlightchannel = sig_channels; % highlight these channels on topoplot
     cfg.highlightcolor = [255/255 192/255 203/255]; % pink colour

     cfg.comment = ['time = [' num2str(cfg.xlim) ']   ' strjoin(stat.label(sig_channels))]; % display time interval & names of sig channels
     %cfg.comment = 'auto'; % display date, xlim (time interval), zlim (amplitude range)
     cfg.commentpos = 'title';   
     %cfg.colorbar = 'yes'; % shows the scaling
     cfg.layout = lay;
     ft_topoplotER(cfg, GA_cue_stayvsswitch);
end  
%}

%% To plot the actual effect (i.e. average ERF of sig channels)
% alt: go to the multiplotER generated earlier, manually select the sig channels & plot

stat = cue_ttype; % can set this to any effect, it's just for reading out the channel labels

% effect in cue window
cfg        = [];
cfg.channel = stat.label(find(cue_ttype.mask)); % autoly retrieve sig channels (only works with cfg.avgovertime = 'yes')
%{'AG017', 'AG018', 'AG019', 'AG022', 'AG023', 'AG025', 'AG029', 'AG063', 'AG064', 'AG143'}; % 10 sig channels in cluster

figure('Name','Average ERF of significant channels - cue window');
ft_singleplotER(cfg, GA_erf.cuechstay, GA_erf.cuechswitch, GA_erf.cueenstay, GA_erf.cueenswitch);
legend(eventnames_real(conds_cue));
xlim([-0.1 0.75]);

% if doing avgovertime, plot vertical lines to indicate the time window
if exist('latency_cue', 'var')
    line([latency_cue(1) latency_cue(1)], ylim, 'Color','black'); % plot a vertical line at start_time
    line([latency_cue(end) latency_cue(end)], ylim, 'Color','black'); % plot a vertical line at end_time
end


% effect in target window
cfg        = [];
cfg.channel = stat.label(find(target_ttype.mask)); % autoly retrieve sig channels (only works with cfg.avgovertime = 'yes')

figure('Name','Average ERF of significant channels - target window');
ft_singleplotER(cfg, GA_erf.targetchstay, GA_erf.targetchswitch, GA_erf.targetenstay, GA_erf.targetenswitch);
legend(eventnames_real(conds_target));
xlim([-0.1 0.75]);

% if doing avgovertime, plot vertical lines to indicate the time window
if exist('latency_target', 'var')
    line([latency_target(1) latency_target(1)], ylim, 'Color','black'); % plot a vertical line at start_time
    line([latency_target(end) latency_target(end)], ylim, 'Color','black'); % plot a vertical line at end_time
end


% Ref code: copied directly from FT_compare_conditions.m
% see also: http://www.fieldtriptoolbox.org/tutorial/cluster_permutation_timelock#the_format_of_the_output
%{

neg_cluster_pvals=[];
pos_cluster_pvals=[];

figure;% % plot negative
cfg = [];
cfg.comment = 'no';
cfg.layout = layout;

cfg.xlim=latency;
cfg.highlight = 'off';

subplot(2,4,1)
eval(['ft_topoplotER(cfg, GM_meg_',cond1,');']);
subplot(2,4,2)
eval(['ft_topoplotER(cfg, GM_meg_',cond2,');']);
subplot(2,4,5)
eval(['ft_topoplotER(cfg, GM_meg_',cond1,'_planar);']);
subplot(2,4,6)
eval(['ft_topoplotER(cfg, GM_meg_',cond2,'_planar);']);

if isfield(eval([cond1,'_vs_',cond2,'_stat']),'posclusters') && ~isempty(eval([cond1,'_vs_',cond2,'_stat.posclusters']))
    eval(['pos_cluster_pvals = [',cond1,'_vs_',cond2,'_stat.posclusters(:).prob];']);
    pos_signif_clust = find(pos_cluster_pvals < 0.05);
    eval(['pos = ismember(',cond1,'_vs_',cond2,'_stat.posclusterslabelmat, pos_signif_clust);']);
    eval(['poscluster_p=([',cond1,'_vs_',cond2,'_stat.posclusters.prob])']);
    cfg.highlight = 'on';
    cfg.highlightchannel = find(pos);
    subplot(3,4,4)
    eval(['plot(GM_meg_',cond1,'.time,mean(GM_meg_',cond1,'.avg(logical(',cond1,'_vs_',cond2,'_stat.posclusterslabelmat),:)))'])
    xlim([-0.2 0.5])
    xlabel('Time (s)')
    ylabel('Amplitude (fT)')
    hold on
    eval(['plot(GM_meg_',cond2,'.time,mean(GM_meg_',cond2,'.avg(logical(',cond1,'_vs_',cond2,'_stat.posclusterslabelmat),:)))'])
else
    sprintf('no significant pos clusters')
end
subplot(2,4,3)
eval(['ft_topoplotER(cfg, GM_meg_',cond1,'_vs_GM_meg_',cond2,');']);
cfg.highlight = 'off';
subplot(2,4,7)
eval(['ft_topoplotER(cfg, GM_meg_planar_',cond1,'_vs_GM_meg_planar_',cond2,');']);
set(gcf, 'Color', 'w');
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure. 

figure;% % plot positive
cfg = [];
cfg.comment = 'no';
cfg.layout = layout;

cfg.xlim=latency;
cfg.highlight = 'off';

subplot(2,4,1)
eval(['ft_topoplotER(cfg, GM_meg_',cond1,');']);
subplot(2,4,2)
eval(['ft_topoplotER(cfg, GM_meg_',cond2,');']);
subplot(2,4,5)
eval(['ft_topoplotER(cfg, GM_meg_',cond1,'_planar);']);
subplot(2,4,6)
eval(['ft_topoplotER(cfg, GM_meg_',cond2,'_planar);']);

if isfield(eval([cond1,'_vs_',cond2,'_stat']),'negclusters') && ~isempty(eval([cond1,'_vs_',cond2,'_stat.negclusters']))
    eval(['neg_cluster_pvals = [',cond1,'_vs_',cond2,'_stat.negclusters(:).prob];']);
    neg_signif_clust = find(neg_cluster_pvals < 0.05);
    eval(['neg = ismember(',cond1,'_vs_',cond2,'_stat.negclusterslabelmat, neg_signif_clust);']);
    eval(['negcluster_p=([',cond1,'_vs_',cond2,'_stat.negclusters.prob])']);
    cfg.highlight = 'on';
    cfg.highlightchannel = find(neg);
    subplot(3,4,4)
    eval(['plot(GM_meg_',cond1,'.time,mean(GM_meg_',cond1,'.avg(logical(',cond1,'_vs_',cond2,'_stat.negclusterslabelmat),:)))'])
    xlim([-0.2 0.5])
    xlabel('Time (s)')
    ylabel('Amplitude (fT)')
    hold on
    eval(['plot(GM_meg_',cond2,'.time,mean(GM_meg_',cond2,'.avg(logical(',cond1,'_vs_',cond2,'_stat.negclusterslabelmat),:)))'])
else
    sprintf('no significant neg clusters')
end
subplot(2,4,3)
eval(['ft_topoplotER(cfg, GM_meg_',cond1,'_vs_GM_meg_',cond2,');']);
cfg.highlight = 'off';
subplot(2,4,7)
eval(['ft_topoplotER(cfg, GM_meg_planar_',cond1,'_vs_GM_meg_planar_',cond2,');']);
set(gcf, 'Color', 'w');
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure. 
end

%}

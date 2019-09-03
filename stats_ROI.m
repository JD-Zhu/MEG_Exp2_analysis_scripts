%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% stats_ROI.m
%
% Grand average & statistical analysis on reconstructed ROI activities
%
%    Q1: do stats across subjects?
%    A: yes. Plot grand ave first, just to get an idea of what effect is there. Then run stats.
%
%    Q2. use ft_timelockstatistics (exactly the same way as doing stats on erf)?
%    A: yes.
%    (but ft_timelockstats takes in erf structures, here we don't have the .avg field, just naked data)
%    A: either hack it into that format, or use EEGlab (statcond, std_stat).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;

% run the #define section
global conds_cue; global conds_target; global eventnames_8;
global ResultsFolder_ROI; % all subjects' ROI data are stored here
global PLOT_SHADE; global colours; % for plotting shaded boundary on each time course
common();


%% Read data

% find all .mat files in ResultsFolder_ROI
files = dir([ResultsFolder_ROI '*_ROI.mat']);

% each cycle reads in one '.mat' file (ie. one subject's ROI results)
for i = 1:length(files)
    filename = [ResultsFolder_ROI files(i).name];
    load(filename);
    allSubjects_ROIs_bySubjects(i) = ROI_activity;
end

% get a list of all the ROI labels
ROIs_label = fieldnames(allSubjects_ROIs_bySubjects(1));

% reformat allSubjects_ROIs: Subject|ROI|condition -> ROI|condition|Subjects
for k = 1:length(ROIs_label)
    ROI_name = ROIs_label{k};
    allSubjects_ROIs.(ROI_name) = allSubjects_reformat(allSubjects_ROIs_bySubjects, ROI_name, eventnames_8);
%{
    for j = 1:length(eventnames_8) % 4 conditions in cue & 4 conditions in target (total 8)
       % if exist('allSubjects_ROI') % if var already exists, append to it
        %if ~isempty(['allSubjects_ROI.' ROIs_label{k} '.' (eventnames_8{j})]) % if var already exists, append to it
        if isfield(allSubjects_ROI, ROIs_label{k}) % if var already exists, append to it
            allSubjects_ROI.(ROIs_label{k}).(eventnames_8{j}) = [allSubjects_ROI.(ROIs_label{k}).(eventnames_8{j}) ROI_activity.(ROIs_label{k}).(eventnames_8{j})];
        else % if first time, simply assign to the first item
            allSubjects_ROI.(ROIs_label{k}).(eventnames_8{j}) = ROI_activity.(ROIs_label{k}).(eventnames_8{j});
        end
    end
%}
end


%% Grand average across subjects

fprintf('\n= COMPUTING & PLOTTING CROSS-SUBJECT AVERAGES =\n');

% each cycle processes one ROI
for k = 1:length(ROIs_label)
    ROI_name = ROIs_label{k};
    fprintf(['\nROI: ' ROI_name '\n']);

    % compute grand average (across subjects) for each condition
    cfg = [];
    cfg.channel   = {'all'}; % there is only one channel (i.e. the virtual sensor for this ROI)
    cfg.latency   = 'all';
    cfg.parameter = 'avg';
    for j = 1:length(eventnames_8)
        cfg.keepindividual = 'no'; % average across subjects
        GA_avg.(eventnames_8{j}) = ft_timelockgrandaverage(cfg, allSubjects_ROIs.(ROI_name).(eventnames_8{j}){:}); 
        
        cfg.keepindividual = 'yes'; % do not average across subjects, keep the data for each individual subject
        GA_keepindi.(eventnames_8{j}) = ft_timelockgrandaverage(cfg, allSubjects_ROIs.(ROI_name).(eventnames_8{j}){:});

        % "{:}" means to use data from all elements of the variable
    end

    GA.(ROI_name) = GA_avg; % store it in the correct field
    GA_indi.(ROI_name) = GA_keepindi; % store it in the correct field
    
    % Plot the GAs
    %{
    % cue-locked
    figure('Name', ['Cue window: GA in ' ROI_name]); hold on
    for j = conds_cue
        plot(GA_erf.(eventnames_8{j}).time, GA_erf.(eventnames_8{j}).avg);
    end
    xlim([-0.1 0.75]);
    legend(eventnames_8(conds_cue));

    % target-locked 
    figure('Name', ['Target window: GA in ' ROI_name]); hold on
    for j = conds_target
        plot(GA_erf.(eventnames_8{j}).time, GA_erf.(eventnames_8{j}).avg);
    end
    xlim([-0.1 0.75]);
    legend(eventnames_8(conds_target));
    %}
end

save([ResultsFolder_ROI 'GA.mat'], 'GA');
save([ResultsFolder_ROI 'GA_individuals.mat'], 'GA_indi');


%% Statistical analysis (to identify time interval of each effect, i.e. temporal clusters)

fprintf('\n= STATS: CLUSTER-BASED PERMUTATION TESTS =\n');

% each cycle processes one ROI
for k = 1:length(ROIs_label)
    
    ROI_name = ROIs_label{k};
    fprintf(['\nROI: ' ROI_name '\n']);

    data = allSubjects_ROIs.(ROI_name); % data for the current ROI
    
    % set some config for the statistical test
    cfg = [];
    cfg.channel = {'all'}; % there is only one channel (i.e. the virtual sensor for this ROI)
    cfg.avgoverchan = 'yes'; % this is necessary (or else FT will ask for cfg.neighbours)
    
    cfg.latency = [-0.1 0.75]; % time interval over which the experimental 
                         % conditions must be compared (in seconds)
    cfg.avgovertime = 'no'; % if yes, this will average over the entire time window chosen in cfg.latency 
                            % (useful when you want to look at a particular component, e.g. to look at M100,
                            % cfg.latency = [0.08 0.12]; cfg.avgovertime = 'yes'; )

    %load([ResultsFolder_ROI 'neighbours.mat']); % this is the sensor layout - it's the same for all subjects (even same across experiments). So just prepare once & save, then load here
    %cfg.neighbours = neighbours;  % same as defined for the between-trials experiment

    cfg.method = 'montecarlo'; %'analytic';
    cfg.statistic = 'depsamplesT'; %cfg.statistic = 'ft_statfun_indepsamplesT'; OR 'ft_statfun_depsamplesFmultivariate';
    cfg.correctm = 'cluster'; %'no';
    cfg.clusteralpha = 0.05;
    cfg.clusterstatistic = 'maxsum';
    %cfg.clusterstatistic = 'wcm'; cfg.wcm_weight = 1;    

    %cfg.minnbchan = 3; % minimum number of neighbourhood channels required to be significant 
                       % in order to form a cluster 
                       % (default: 0, ie. each single channel can be considered a cluster).
                       % 4 or 5 is a good choice; 2 is too few coz it's even below
                       % the resolution of the sensor layout(??)

    cfg.tail = 0;
    cfg.clustertail = 0; % 2 tailed test
    cfg.alpha = 0.1; % report all effects with p < 0.1
    cfg.correcttail = 'prob'; % correct for 2-tailedness
    cfg.numrandomization = 2000; % Rule of thumb: use 500, and double this number if it turns out 
        % that the p-value differs from the critical alpha-level (0.05 or 0.01) by less than 0.02

    numSubjects = length(files);
    within_design_2x2 = zeros(2, 2*numSubjects);
    within_design_2x2(1, :) = repmat(1:numSubjects, 1, 2);
    within_design_2x2(2, 1:numSubjects) = 1;
    within_design_2x2(2, numSubjects+1:2*numSubjects) = 2;

    cfg.design = within_design_2x2;
    cfg.uvar  = 1; % row of design matrix that contains unit variable (in this case: subjects)
    cfg.ivar  = 2; % row of design matrix that contains independent variable (i.e. the conditions)

    % Run the statistical tests
    
    % Interaction (i.e. calc sw$ in each lang, then submit the 2 sw$ for comparison)
    fprintf('\nCUE window -> Testing lang x ttype interaction:\n');
    [timelock1, timelock2] = combine_conds_for_T_test('fieldtrip', 'interaction', data.cuechstay, data.cuechswitch, data.cueenstay, data.cueenswitch);
    %cfg.latency = latency_cue; % time interval over which the experimental 
    [cue_interaction.(ROI_name)] = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:});
    fprintf('\nTARGET window -> Testing lang x ttype interaction:\n');
    [timelock1, timelock2] = combine_conds_for_T_test('fieldtrip', 'interaction', data.targetchstay, data.targetchswitch, data.targetenstay, data.targetenswitch); %'2-1 vs 4-3');
    %cfg.latency = latency_target; % time interval over which the experimental 
    [target_interaction.(ROI_name)] = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:});
    
    % Main effect of lang (collapse across stay-switch)
    fprintf('\nCUE window -> Main effect of lang:\n');
    [timelock1, timelock2] = combine_conds_for_T_test('fieldtrip', 'main_12vs34', data.cuechstay, data.cuechswitch, data.cueenstay, data.cueenswitch);
    %cfg.latency = latency_cue; % time interval over which the experimental 
    [cue_lang.(ROI_name)] = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:});
    fprintf('\nTARGET window -> Main effect of lang:\n');
    [timelock1, timelock2] = combine_conds_for_T_test('fieldtrip', 'main_12vs34', data.targetchstay, data.targetchswitch, data.targetenstay, data.targetenswitch); %'2-1 vs 4-3');
    %cfg.latency = latency_target; % time interval over which the experimental 
    [target_lang.(ROI_name)] = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:});

    % Main effect of switch (collapse across langs)
    fprintf('\nCUE window -> Main effect of ttype:\n');
    [timelock1, timelock2] = combine_conds_for_T_test('fieldtrip', 'main_13vs24', data.cuechstay, data.cuechswitch, data.cueenstay, data.cueenswitch);
    %cfg.latency = latency_cue; % time interval over which the experimental 
    [cue_ttype.(ROI_name)] = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:});
    fprintf('\nTARGET window -> Main effect of ttype:\n');
    [timelock1, timelock2] = combine_conds_for_T_test('fieldtrip', 'main_13vs24', data.targetchstay, data.targetchswitch, data.targetenstay, data.targetenswitch); %'2-1 vs 4-3');
    %cfg.latency = latency_target; % time interval over which the experimental 
    [target_ttype.(ROI_name)] = ft_timelockstatistics(cfg, timelock1{:}, timelock2{:});

end

save([ResultsFolder_ROI 'stats.mat'], 'cue_interaction', 'cue_lang', 'cue_ttype', 'target_interaction', 'target_lang', 'target_ttype');


%% Find the effects & plot them
% Automatically check all the stats output & read out the time interval
% of each effect (from the stat.mask field)

stats = load([ResultsFolder_ROI 'stats.mat']);
load([ResultsFolder_ROI 'GA.mat']);
load([ResultsFolder_ROI 'GA_individuals.mat']);

% make directory to store the output figures
mkdir([ResultsFolder_ROI 'Figures\\non-sig\\']);

ROIs_names = fieldnames(GA); % get the list of ROI names

% baseline correction b4 plotting
% Note: This is the right place to do baseline correction. We decided not to do it
% b4 running stats (i.e. at single-subject level) - see my email for expla
cfg = [];
cfg.feedback = 'no';
cfg.baseline = [-0.1 0];
for k = 1:length(ROIs_names) % each cycle handles one ROI
    ROI_name = ROIs_names{k};
    for j = 1:length(eventnames_8)
        GA.(ROI_name).(eventnames_8{j}) = ft_timelockbaseline(cfg, GA.(ROI_name).(eventnames_8{j})); 
    end
end

% loop thru all 6 stats output (cue/target lang/ttype/interxn) and loop thru all ROIs in each,
% check whether the .mask field has any non-zero entries 
% (these are the effects & they are already cluster-corrected, so doesn't need to be consecutive 1s)
fprintf('\nThe following effects were detected:\n');
stats_names = fieldnames(stats);
for i = 1:length(stats_names) % each cycle handles one effect (e.g. cue_lang)
    stat_name = stats_names{i};
    ROIs_names = fieldnames(stats.(stat_name)); % get the list of ROI names
    
    for k = 1:length(ROIs_names) % each cycle handles one ROI
        ROI_name = ROIs_names{k};
                        
        % if the .mask contains any 1s, that's an effect (the rest are 0s)
        mask = stats.(stat_name).(ROI_name).mask; 
        effect = find(stats.(stat_name).(ROI_name).mask); 

        % do a GA plot for all contrasts that have an effect (will add gray box to show effect interval later)
        % also do a GA plot for all ROIs that don't have an effect (save in 'non-sig' folder); we don't need to plot the 3 contrasts for each ROI ('lang', 'ttype' & 'interaction') - they are the same, just plot one
        if ( ~isempty(effect) || strcmp(stat_name(end-3:end), 'tion') ) % can't compare the whole word 'interaction' here, coz some stat_names (e.g. 'cue_lang') are a shorter string
            % GA plot
            figure('Name', [stat_name ' in ' ROI_name], 'Position', get(0, 'Screensize')); % make the figure full-screen size
            hold on;
            
            if strcmp(stat_name(1:3), 'cue') % this effect occurs in cue window
                
                % each cycle plots 1 line (ie. 1 condition)
                for j = conds_cue 
                    if strcmp(PLOT_SHADE, 'no') % do not plot shaded boundary, just plot a single line                  
                        plot(GA.(ROI_name).(eventnames_8{j}).time, GA.(ROI_name).(eventnames_8{j}).avg);
                    else % calc the margin for shaded boundary (stdev / sem / CI) at every time point
                        allsubjects = GA_indi.(ROI_name).(eventnames_8{j}).individual;
                        margin = calc_margin(allsubjects, PLOT_SHADE);
                        
                        % plot time course with shaded boundary
                        boundedline(GA.(ROI_name).(eventnames_8{j}).time, GA.(ROI_name).(eventnames_8{j}).avg, margin(:), 'alpha', 'transparency',0.15, colours(j));                        
                    end
                end
                
                xlim([-0.1 0.75]); 
                
            elseif strcmp(stat_name(1:6), 'target') % this effect occurs in target window
                
                for j = conds_target
                    if strcmp(PLOT_SHADE, 'no') % do not plot shaded boundary, just plot a single line                  
                        plot(GA.(ROI_name).(eventnames_8{j}).time, GA.(ROI_name).(eventnames_8{j}).avg);
                    else % calc the margin for shaded boundary (stdev / sem / CI) at every time point
                        allsubjects = GA_indi.(ROI_name).(eventnames_8{j}).individual;
                        margin = calc_margin(allsubjects, PLOT_SHADE);

                        % plot time course with shaded boundary
                        boundedline(GA.(ROI_name).(eventnames_8{j}).time, GA.(ROI_name).(eventnames_8{j}).avg, margin(:), 'alpha', 'transparency',0.15, colours(j));
                    end                    
                end
                
                xlim([-0.1 0.55]); 

            else % should never be here
                fprintf('Error: an effect is found, but its not in either cue nor target window.\n');
            end
            
            
            % set properties for axes, lines, and text
            xlabel('Seconds');
            ylabel('Ampere per square metre');
            set(gca, 'LineWidth',1.5, 'FontSize',22); % set axes properties
            box on; % draw a border around the figure

            % specify the legend manually (otherwise it will include
            % each shaded patch as an item too). For some reason,
            % the order of the lines are reversed when you grab them
            lines = findall(gcf, 'Type','line');
            legend([lines(4) lines(3) lines(2) lines(1)], {'Mandarin (L1) stay', 'Mandarin (L1) switch', 'English (L2) stay', 'English (L2) switch'}, 'Location','northwest', 'FontSize',30);
            set(lines, 'Linewidth',3); % line thickness
                
            % reference code:
            %{
            p = findobj(gcf); % get the handles associated with the current figure
            
            allaxes = findall(p,'Type','axes');
            alllines = findall(p,'Type','line');
            alltext = findall(p,'Type','text');
            
            set(allaxes,'FontName','Arial','FontWeight','Bold','LineWidth',2,'FontSize',14);
            set(alllines,'Linewidth',3);
            set(alltext,'FontName','Arial','FontWeight','Bold','FontSize',14);
            %}


            % if this GA plot is for a non-sig ROI, then save the figure into the 'non-sig' folder
            if isempty(effect)
                
                % for these plots, we want to delete the legends
                hl = findobj(gcf, 'type','legend');
                delete(hl);
                
                % maximise the figure before saving
                %set(gcf, 'Position', get(0, 'Screensize'));

                filename = [ROI_name '_' stat_name(1:3) '.png'];
                %saveas(gcf, [ResultsFolder_ROI 'Figures\\non-sig\\' filename]); % this fn does not maintain the aspect ratio, font size, etc
                export_fig(gcf, [ResultsFolder_ROI 'Figures\\non-sig\\' filename]); % use this tool to save the figure exactly as shown on screen

            else % if there is any effect present, find all the clusters so we can output to console & mark on the plot

                % find the start & end of each cluster
                start_points = find(diff(mask) == 1); % transitions from 0 to 1
                end_points = find(diff(mask) == -1); % transitions from 1 to 0

                % check if we are missing the very first start point
                if (mask(1) == 1)
                    start_points = [1 start_points];
                end
                % check if we are missing the very last end point
                if (mask(end) == 1)
                    end_points = [end_points length(mask)];
                end
                % sanity check
                assert (length(start_points) == length(end_points));

                % increase the index for all start points by 1 (coz the "transition" found by diff 
                % is the position of the last '0' before it turns into '1')
                start_points = start_points + 1;


                % produce console output for each cluster & mark it on the GA plot
                fprintf('%s has an effect in %s:\n', ROI_name, stat_name);

                for cluster = 1:length(start_points) % start_points contains a list of starting points (one for each cluster)
                    start_sample = start_points(cluster);
                    end_sample = end_points(cluster);

                    % read out the necessary info
                    pvalue = stats.(stat_name).(ROI_name).prob(start_sample); % p-value is the same for all time points in a cluster, so we just read it from the first time point
                    start_time = stats.(stat_name).(ROI_name).time(start_sample);
                    end_time = stats.(stat_name).(ROI_name).time(end_sample); 

                    fprintf('    p = %.4f, between %.f~%.f ms (significant at samples %d to %d).\n', ...
                        pvalue, start_time*1000, end_time*1000, start_sample, end_sample); % convert units to ms

                    % mark the cluster interval on the GA plot
                    %line([start_time start_time], ylim, 'Color','black'); % plot a vertical line at start_time
                    %line([end_time end_time], ylim, 'Color','black'); % plot a vertical line at end_time

                    % create shaded region indicating effect duration
                    ylimits = ylim; ylow = ylimits(1); yhigh = ylimits(2);
                    x = [start_time end_time end_time start_time]; % specify x,y coordinates of the 4 corners
                    y = [ylow ylow yhigh yhigh];
                    patch(x,y,'black', 'FaceAlpha',0.15) % draw the shade (FaceAlpha is transparency)
                    ylim(ylimits); % ensure ylim doesn't get expanded
                end

                % maximise the figure before saving
                %set(gcf, 'Position', get(0, 'Screensize'));
                
                % save the figure
                filename = [ROI_name '_' stat_name '.png'];
                %saveas(gcf, [ResultsFolder_ROI 'Figures\\' filename]); % this fn does not maintain the aspect ratio, font size, etc
                export_fig(gcf,[ResultsFolder_ROI 'Figures\\' filename]); % use this tool to save the figure exactly as shown on screen
                
                % old code to check multiple clusters
                %{
                % if the .mask contains any non-zero entries, that's an effect
                effect = find(stats.(stat_name).(ROI_name).mask); 
                if ~isempty(effect) % if there is an effect, we print it out & plot the ROI timecourse for each cond

                % check whether there are any jumps in the indices (i.e. multiple clusters)
                jump_positions = find(diff(effect) ~= 1);

                % find the start point for the first temporal cluster
                start_sample = effect(1);

                % find the end point for each cluster & print out this cluster
                for cluster = 1:length(jump_positions)

                    end_sample = effect(jump_positions(cluster));

                    % read out the necessary info
                    pvalue = stats.(stat_name).(ROI_name).prob(start_sample); % p-value is the same for all time points in a cluster, so we just read it from the first time point
                    start_time = stats.(stat_name).(ROI_name).time(start_sample);
                    end_time = stats.(stat_name).(ROI_name).time(end_sample); 

                    fprintf('%s has an effect in %s (p = %.4f), between %.f~%.f ms (significant at samples %d to %d).\n', ROI_name, stat_name, pvalue, start_time*1000, end_time*1000, start_sample, end_sample); % convert units to ms

                    % start point for the next cluster
                    start_sample = effect(jump_positions(cluster) + 1);
                end

                % find the end point for the final cluster & print out this cluster
                end_sample = effect(end); 
                % read out the necessary info
                pvalue = stats.(stat_name).(ROI_name).prob(start_sample); % p-value is the same for all time points in a cluster, so we just read it from the first time point
                start_time = stats.(stat_name).(ROI_name).time(start_sample);
                end_time = stats.(stat_name).(ROI_name).time(end_sample); 

                fprintf('%s has an effect in %s (p = %.4f), between %.f~%.f ms (significant at samples %d to %d).\n', ROI_name, stat_name, pvalue, start_time*1000, end_time*1000, start_sample, end_sample); % convert units to ms

                %}
            end
            
            hold off;
        end
    end
end

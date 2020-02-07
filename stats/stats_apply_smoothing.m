% Script to apply smoothing on all single-subject ROI timecourses, 
% to make them less "spiky" (could help to detect effects more easily)
%
% Author: Judy Zhu (github.com/JD-Zhu)
%
% 
%% setup
global ResultsFolder_ROI; % all subjects' ROI data are stored here
common();

% SPECIFY location containing the original (ie. not yet smoothed) timecourses
original_folder = [ResultsFolder_ROI 'TSPCA10000_3_freeori'];


%% start
% create a folder to store the smoothed timecourses
smoothed_folder = [original_folder '_smooth\\'];
mkdir(smoothed_folder);
original_folder = [original_folder '\\'];

% find all .mat files in original_folder
files = dir([original_folder '*_ROI.mat']);

% each cycle reads in one '.mat' file (ie. one subject's ROI timecourses)
for i = 1:length(files)
    filename = [original_folder files(i).name];
    load(filename);
    
    % loop thru each ROI
    ROIs_label = fieldnames(ROI_activity);
    for k = 1:length(ROIs_label)
        ROI_name = ROIs_label{k};
        
        % loop thru each cond
        conds_label = fieldnames(ROI_activity.(ROI_name));
        for j = 1:length(conds_label)
            cond_name = conds_label{j};
            
            % read the original (ie. not yet smoothed) timecourse
            timecourse = ROI_activity.(ROI_name).(cond_name).avg;
            
            % smooth it, then replace the original
            timecourse_smooth = smooth(timecourse);
            ROI_activity.(ROI_name).(cond_name).avg = timecourse_smooth';
            
            % plot for sanity check
            %{
            figure;
            plot(ROI_activity.(ROI_name).(cond_name).time, timecourse);
            hold on;
            plot(ROI_activity.(ROI_name).(cond_name).time, timecourse_smooth);
            hold off;
            %}
        end
    end
    
    % save the smoothed timecourses to new location
    save([smoothed_folder files(i).name], 'ROI_activity');
end

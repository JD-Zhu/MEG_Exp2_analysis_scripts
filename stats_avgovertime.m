% First, run the top 2 sections in stats_ROI.m



% Then, run the following code to extract a single value for each 
% subject in each condition (avg over time)

% Based on the finding:
% Bi_sw in LIFG, p = 0.0770, between 425~450 ms
ROI_name = 'LIFG';
start_sample = 246; % check the .time field in "data", to find the sample indices corresponding to the time interval above
end_sample = 251;

data = allSubjects_ROIs.(ROI_name); % data for the current ROI
T = table(); % create an empty table (to store the output)
    
cond_names = fieldnames(data);
for j = 1:length(cond_names) % each cycle handles one cond (e.g. NatStay)
    cond_name = cond_names{j};
    allsubjects = data.(cond_name); % all subjects' time courses for this cond
    
    for i = 1:length(allsubjects) % each cycle handles one subject
        subject_timecourse = allsubjects{i}.avg;
        avg = mean(subject_timecourse(start_sample:end_sample)); % take the average over the pre-defined time interval
        T.(cond_name)(i,1) = avg; % store into the appropriate cell in the table
    end   
end

% save to csv file
output_file = [ResultsFolder_ROI_thisrun 'avgovertime\\' ROI_name '_' start_sample '-' end_sample 'ms.csv'];
writetable(T, output_file);

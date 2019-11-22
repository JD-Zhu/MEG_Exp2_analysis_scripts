% Which model? (Full factorial / Flexible factorial)
which_model = 'flexible'; % 'full';

% Select input folder (where to read in the images)
image_dir = 'E:\Judy\Exp2\6_MEG-data\Results_ROI\TSPCA10000_3\SPM_temp\';

% Select output folder (output files will be saved here)
results_dir = 'E:\Judy\Exp2\6_MEG-data\Results_ROI\TSPCA10000_3\SPM_results\';

%=========================================================================%

load('ROIs_label.mat');

% Each cycle processes one ROI
for k = 1:length(ROIs_label)
    ROI_name = ROIs_label{k};
    
    % Generate list of images (required input for "Specify 2nd-level")
    if strcmp(which_model, 'flexible')
        % generate images for all subjects in all conds (for this ROI)
        conds = genImageList_FlexibleFactorial(image_dir, ROI_name); 
        save([results_dir 'conds.mat'], 'conds');
        
        % generate the design matrix
        design_complete = genDesignMatrix();
        save([results_dir 'design.mat'], 'design_complete');

    elseif strcmp(which_model, 'full')
        % generate list of 24 images (one image per subject) in each 3x3 condition
        cond_1_1 = genImageList_FullFactorial(image_dir, ROI_name, 'NatSingle'); % 24 images in the NatSingle cond
        cond_1_2 = genImageList_FullFactorial(image_dir, ROI_name, 'NatStay');
        cond_1_3 = genImageList_FullFactorial(image_dir, ROI_name, 'NatSwitch');
        cond_2_1 = genImageList_FullFactorial(image_dir, ROI_name, 'ArtSingle');
        cond_2_2 = genImageList_FullFactorial(image_dir, ROI_name, 'ArtStay');
        cond_2_3 = genImageList_FullFactorial(image_dir, ROI_name, 'ArtSwitch');
        cond_3_1 = genImageList_FullFactorial(image_dir, ROI_name, 'BiSingle');
        cond_3_2 = genImageList_FullFactorial(image_dir, ROI_name, 'BiStay');
        cond_3_3 = genImageList_FullFactorial(image_dir, ROI_name, 'BiSwitch');

        save([results_dir 'conds.mat'], 'cond_1_1', 'cond_1_2', 'cond_1_3', 'cond_2_1', 'cond_2_2', 'cond_2_3', 'cond_3_1', 'cond_3_2', 'cond_3_3');
    end


    % = The 2 steps below are based on the GUI steps = %
    % ("Specify 2nd-level" & "Estimate")
    
    % Specify 2nd-level model
    nrun = 1; % enter the number of runs here
    if strcmp(which_model, 'flexible')
        jobfile = {'SPM_specify_model_job.m'};
    elseif strcmp(which_model, 'full')
        jobfile = {'SPM_specify_model_FullFactorial_job.m'};
    end
    jobs = repmat(jobfile, 1, nrun);
    inputs = cell(0, nrun);
    for crun = 1:nrun
    end
    spm('defaults', 'EEG');
    spm_jobman('run', jobs, inputs{:});
    
    % Estimate model
    nrun = 1; % enter the number of runs here
    jobfile = {'SPM_estimate_model_job.m'};
    jobs = repmat(jobfile, 1, nrun);
    inputs = cell(0, nrun);
    for crun = 1:nrun
    end
    spm('defaults', 'EEG');
    spm_jobman('run', jobs, inputs{:});
    
    
    % move all model results for this ROI into a separate folder
    files = dir(results_dir);  %get all files/folders
    isfile= ~[files.isdir]; %determine index of files vs folders
    filenames = {files(isfile).name}; %create cell array of file names
    curr_ROI_folder = [results_dir ROI_name];
    mkdir(curr_ROI_folder);
    for i = 1:length(filenames)
        movefile([results_dir filenames{i}], curr_ROI_folder);
    end
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                       Sub-functions                             %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate list of images for all subjects in all conds (for the given ROI)
function images = genImageList_FlexibleFactorial(image_dir, ROI_name)
    images = {};
    NumSubjects = 24;
    cond_names = {'NatSingle', 'NatStay', 'NatSwitch', 'ArtSingle', 'ArtStay', 'ArtSwitch', 'BiSingle', 'BiStay', 'BiSwitch'};
    
    for i = 1:NumSubjects
        for j = 1:length(cond_names)
            image_file = [image_dir ROI_name '_' cond_names{j} '_subj' int2str(i) '_evoked\condition_Undefined.nii,1'];
            images{end+1} = image_file; % append this image_file to the end of the list
        end
    end
    
    images = images';
end

% generate list of images for all subjects in the given cond_name (for the given ROI)
function images = genImageList_FullFactorial(image_dir, ROI_name, cond_name)
    images = [];
    NumSubjects = 24;
    
    for i = 1:NumSubjects
        image_file = [image_dir ROI_name '_' cond_name '_subj' int2str(i) '_evoked\condition_Undefined.nii,1'];
        images{i} = image_file;
    end
end

% generate the full design matrix for flexible-factorial model
function design_complete = genDesignMatrix()
    design_complete = [];

    design =   [1 1
                1 2
                1 3
                2 1
                2 2
                2 3
                3 1
                3 2
                3 3];

    NumConds = size(design, 1);
    NumSubjects = 24;

    % each cycle creates the 9 rows for one subject
    for i = 1:NumSubjects
        %col_1 = ones(NumConds, 1); % column 1 is all 1s
        col_1 = repmat(1, [NumConds 1]); % column 1 is all 1s
        col_2 = repmat(i, [NumConds 1]); % column 2 is the subject number
        design_oneSubject = [col_1 col_2 design];

        design_complete = vertcat(design_complete, design_oneSubject);
    end
end

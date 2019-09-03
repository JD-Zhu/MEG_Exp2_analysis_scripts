%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MRI Estimation for MEG Sourcespace (MEMES)
%
% find a most appropriate structural MRI based on the individual's headshape
% and then coregister it with the MEG data and headshape information
%
% INPUTS: (via param passing)
% - dir_name        = directory name for the output of your coreg
% - coreg_output    = where to store the output from MEMES
% - confile         = full path to the con file
% - mrkfile         = full path to the mrk file
% - mri_file        = full path to the NIFTI structural MRI file
% - hspfile         = full path to the hsp (polhemus headshape) file
% - elpfile         = full path to the elp file
% - hsp_points      = number of points for downsampling the headshape (try
% 100-200) -> option no longer avail, specified by magic number
% - scalpthreshold  = threshold for scalp extraction (try 0.05 if unsure) -> option no longer avail, specified by magic number
% - bad_coil        = list of bad coils (up to length of 2). Enter as:
%                         {LPAred','RPAyel','PFblue','LPFwh','RPFblack'}
%
% OUTPUTS: (via save)
% - grad_trans              = correctly aligned sensor layout
% - headshape_downsampled   = downsampled headshape (original variable name I know)
% - mri_realigned           = the mri realigned based on fiducial points
% - trans_matrix            = transformation matrix for accurate coregistration
% - headmodel_singleshell   = coregistered singleshell headmodel
%
% Written by Robert Seymour Oct 2017 (some subfunctions written by Paul
% Sowman)
%
% Original script:
% https://github.com/Macquarie-MEG-Research/MEMES
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MEMES(dir_name,coreg_output,elpfile,hspfile,confile,mrkfile,path_to_MRI_library,mesh_library,initial_mri_realign,bad_coil)

    fprintf('\nThis is MEMES.m\n');
    
    %% Check inputs
    disp('Performing input check');
    assert(length(bad_coil)<3,'You need at least 3 good coils for accurate alignment. Also make sure you enter bad_coil strings in curly brackets {}');
    % If Path to MRI library doesn't end with / or \ throw up and error
    if ismember(path_to_MRI_library(end),['/','\']) == 0
        error('!!! Path to MRI library must end with / or \ !!!');
    end

    % Check if bad_coils are entered correctly
    if strcmp(bad_coil,'')
        disp('No bad coils marked');
    else
        for check1 = 1:length(bad_coil)
            if ismember(bad_coil{check1},{'','LPAred','RPAyel','PFblue','LPFwh','RPFblack'}) == 0
                error('!!! Please enter bad_coils correctly in the form {LPAred,RPAyel,PFblue,LPFwh,RPFblack} !!!');
            end
        end
    end

    
    %% list of HCP subjects
    subject = {'100307';'102816';'104012';'105923';'106521';'108323';...
        '109123';'111514';'112920';'113922';'116524';'116726';'125525';...
        '133019';'140117';'146129';'149741';'151526';'153732';'154532';...
        '156334';'158136';'162026';'162935';'164636';'166438';'169040';...
        '172029';'174841';'175237';'175540';'177746';'179245';'181232';...
        '182840';'185442';'187547';'189349';'191033';'191437';'191841';...
        '192641';'195041';'198653';'200109';'204521';'205119';'212318';...
        '212823';'214524';'221319';'223929';'233326';'248339';'250427';...
        '255639';'257845';'283543';'287248';'293748';'352132';'352738';...
        '353740';'358144';'406836';'433839';'500222';'512835';'555348';...
        '559053';'568963';'581450';'599671';'601127';'660951';'662551';...
        '665254';'667056';'679770';'680957';'706040';'707749';'715950';...
        '725751';'735148';'783462';'814649';'825048';'872764';'877168';...
        '891667';'898176';'912447';'917255';'990366'};

    nii_filename = '\\MEG\\anatomy\\T1w_acpc_dc_restore.nii'; % same for all subjects

    % make folder to contain the MEMES output
    if ~exist(coreg_output)
        mkdir(coreg_output);
    end

    % CD to right place
    cd(dir_name); fprintf('\n CDd to the right place\n');

    %addpath('/Users/44737483/Documents/scripts_mcq/alien');

    % Get Polhemus Points
    [shape] = parsePolhemus(elpfile,hspfile);

    % Read the grads from the con file
    grad_con                    = ft_read_sens(confile); %in cm, load grads

    % Read mrk_file
    mrk      = ft_read_headshape(mrkfile,'format','yokogawa_mrk');
    mrk      = ft_convert_units(mrk,'cm'); %in cm
  
    %% Perform Realignment Using Paul's Fancy Functions
    if strcmp(bad_coil,'')
        disp('NO BAD MARKERS');
        markers                     = mrk.fid.pos([2 3 1 4 5],:);%reorder mrk to match order in shape
        [R,T,Yf,Err]                = rot3dfit(markers,shape.fid.pnt(4:end,:));%calc rotation transform
        meg2head_transm             = [[R;T]'; 0 0 0 1];%reorganise and make 4*4 transformation matrix

        disp('Performing re-alignment');
        grad_trans                  = ft_transform_geometry_PFS_hacked(meg2head_transm,grad_con); %Use my hacked version of the ft function - accuracy checking removed not sure if this is good or not
        grad_trans.fid              = shape; %add in the head information
        save ([coreg_output 'grad_trans.mat'], 'grad_trans');

    % Else if there is a bad marker
    else
        fprintf(''); disp('TAKING OUT BAD MARKER(S)');

        badcoilpos = [];

        % Identify the bad coil
        for num_bad_coil = 1:length(bad_coil)
            pos_of_bad_coil = find(ismember(shape.fid.label,bad_coil{num_bad_coil}))-3;
            badcoilpos(num_bad_coil) = pos_of_bad_coil;
        end

        % Re-order mrk file to match elp file
        markers               = mrk.fid.pos([2 3 1 4 5],:);%reorder mrk to match order in shape
        % Now take out the bad marker(s) when you realign
        markers(badcoilpos,:) = [];

        % Get marker positions from elp file
        fids_2_use = shape.fid.pnt(4:end,:);
        % Now take out the bad marker(s) when you realign
        fids_2_use(badcoilpos,:) = [];

        % If there are two bad coils use the ICP method, if only one use
        % rot3dfit as usual
        disp('Performing re-alignment');

        if length(bad_coil) == 2
            fprintf('\nTWO BAD COILS!!!!!!!!!!\n\n');
            
            [R, T, err, dummy, info]    = icp(fids_2_use', markers','Minimize', 'point');
            meg2head_transm             = [[R T]; 0 0 0 1];%reorganise and make 4*4 transformation matrix
            grad_trans                  = ft_transform_geometry_PFS_hacked(meg2head_transm,grad_con); %Use my hacked version of the ft function - accuracy checking removed not sure if this is good or not
            grad_trans.fid              = shape; %add in the head information
        else
            [R,T,Yf,Err]                = rot3dfit(markers,fids_2_use);%calc rotation transform
            meg2head_transm             = [[R;T]'; 0 0 0 1];%reorganise and make 4*4 transformation matrix
            grad_trans                  = ft_transform_geometry_PFS_hacked(meg2head_transm,grad_con); %Use my hacked version of the ft function - accuracy checking removed not sure if this is good or not
            grad_trans.fid              = shape; %add in the head information
        end
    end

    % Create figure to view relignment
    hfig = figure;
    subplot(2,2,1);ft_plot_headshape(shape);
    hold on; ft_plot_sens(grad_trans); view([180, 0]);
    subplot(2,2,2);ft_plot_headshape(shape);
    hold on; ft_plot_sens(grad_trans); view([-90, 0]);
    subplot(2,2,3);ft_plot_headshape(shape);
    hold on; ft_plot_sens(grad_trans); view([0, 0]);
    hax = subplot(2,2,4);ft_plot_headshape(shape);
    hold on; ft_plot_sens(grad_trans); view([90, 0]);

    
    % Get headshape downsampled to 100 points with facial info preserved
    headshape_downsampled = downsample_headshape_noface(hspfile,100,grad_trans);

    % Rotate sensors and headshape about z-axis
    rot180mat       = rotate_about_z(180);
    grad_trans      = ft_transform_geometry(rot180mat,grad_trans);
    headshape_downsampled = ft_transform_geometry(rot180mat,headshape_downsampled);

    
    %% Perform ICP

    % Initialise coreg error (ORE) to 1, for all candidate scalp surfaces
    error_term = zeros(1, length(mesh_library)); 
    % Variable to hold the transformation matrices
    trans_matrix_library = [];

    % loop thru all scalp surfaces, coregister each to digitised headshape
    % and compute ORE (smallest error = winner)
    for m = 1:length(mesh_library)
        % Perform ICP (interative closest point) to compute ORE & trans matrix
        numiter = 50; 
        [R, t, err, dummy, info] = icp(mesh_library{m}.pos', headshape_downsampled.pos', numiter, 'Minimize', 'plane', 'Extrapolation', true,'WorstRejection', 0.05);

        % Add ORE for this candidate scalp to the list of error values
        error_term(m) = err(end);

        % Add transformation matrix for this candidate to the list
        trans_matrix_library{m} = inv([real(R) real(t); 0 0 0 1]);

        fprintf('Completed iteration %d of %d\n', m, length(mesh_library));
    end

    %% Make pretty figure
    error_term_sorted = sort(error_term, 'descend');
    losers = find(ismember(error_term,error_term_sorted(1:3))); % worst 3 examples
    middles = find(ismember(error_term,error_term_sorted(46:48))); % middle 3 examples
    winners = find(ismember(error_term,error_term_sorted(end-2:end))); % best 3 examples

    concat = [winners middles losers];

    % Create figure to summarise the losers,middles and winners
    figure;
    for i = 1:9

        mesh_spare = mesh_library{(concat(i))};
        mesh_spare.pos = ft_warp_apply(trans_matrix_library{(concat(i))}, mesh_spare.pos);

        subplot(3,3,i)
        ft_plot_mesh(mesh_spare,'facecolor',[238,206,179]./255,'EdgeColor','none','facealpha',0.8); hold on;
        camlight; hold on; view([-180,-10]);
        if ismember(i,1:3)
            title(sprintf('BEST: %d', error_term((concat(i)))));
        elseif ismember(i,4:6)
            title(sprintf('MIDDLE: %d', error_term((concat(i)))));
        elseif ismember(i,7:9)
            title(sprintf('WORST: %d', error_term((concat(i)))));
        end

        ft_plot_headshape(headshape_downsampled);

        if i == 9
            print([coreg_output 'best_middle_worst_examples'],'-dpdf','-r200');
        end
    end

    %% Use the best for to create a source model for MEG source analysis

    winner = find(error_term == min(error_term));
    fprintf('\nThe winning MRI is number %d of %d\n',winner,length(mesh_library));
    trans_matrix = trans_matrix_library{winner};

    % Create figure to show ICP fit
    mesh_spare = mesh_library{winner};
    mesh_spare.pos = ft_warp_apply(trans_matrix, mesh_spare.pos);

    figure;ft_plot_mesh(mesh_spare,'facecolor',[238,206,179]./255,'EdgeColor','none','facealpha',0.8); hold on;
    camlight; hold on; view([-180,-10]);
    title(error_term(winner));
    ft_plot_headshape(headshape_downsampled);

    print([coreg_output 'winning_sourcemodel'],'-dpdf','-r200');

    try
        % % Make fancy video
        c = datestr(clock); %time and date

        figure;
        ft_plot_mesh(mesh_spare,'facecolor',[238,206,179]./255,'EdgeColor','none','facealpha',0.8); hold on;
        camlight; hold on;
        ft_plot_headshape(headshape_downsampled); title(sprintf('%s.   Error of ICP fit = %d' , c, error_term(winner)));
        OptionZ.FrameRate=15;OptionZ.Duration=5.5;OptionZ.Periodic=true;
        CaptureFigVid([0,0; 360,0], 'ICP_quality',OptionZ)

    catch
        fprintf('You need CaptureFigVid in your path for fancy videos\n');
    end

    % Get MRI of winning subject
    mri_file = [path_to_MRI_library subject{winner} nii_filename];

    mri_orig                    = ft_read_mri(mri_file); % in mm, read in mri from DICOM
    mri_orig = ft_convert_units(mri_orig,'cm'); mri_orig.coordsys = 'neuromag';

    mri_orig.transform = initial_mri_realign{winner};
    mri_realigned = mri_orig; clear mri_orig;

    
    % at this stage, save some necessary variables, because ft_volumesegment sometimes
    % crashes (and when it does, the variables are wiped & we have to start
    % again from the beginning of this fn)
    mkdir ([coreg_output 'tmp\\']);
    save ([coreg_output 'tmp\\tmp_save.mat'], 'coreg_output', 'headshape_downsampled', 'trans_matrix', 'grad_trans'); 
    %load ([coreg_output 'tmp\\tmp_save.mat']);
    
    
    %%
    % Segment
    fprintf('\nSegmenting the MRI... This may take a while...\n');
    
    ft_defaults % reset paths, or else ft_volumesegment will crash
    cfg           = [];
    cfg.output    = 'brain';
    try   % use try-catch, just in case it crashes
        mri_segmented  = ft_volumesegment(cfg, mri_realigned);
    catch
        fprintf('\nNote: ft_volumesegment crashed. Giving it one more go...\n');
        
        % if failed, try again
        ft_defaults
        cfg           = [];
        cfg.output    = 'brain';
        mri_segmented  = ft_volumesegment(cfg, mri_realigned);
    end

    % Create singleshell headmodel
    cfg = [];
    cfg.method = 'singleshell';

    headmodel_singleshell = ft_prepare_headmodel(cfg, mri_segmented); % in cm, create headmodel

    % Apply transformation matrix
    headmodel_singleshell.bnd.pos = ft_warp_apply(trans_matrix, headmodel_singleshell.bnd.pos);
    mri_realigned = ft_transform_geometry(trans_matrix, mri_realigned);
    
    % sanity check (plot mri & headmodel together)
    %ft_determine_coordsys(mri, 'interactive','no'); hold on
    %ft_plot_vol(headmodel);


    figure;ft_plot_headshape(headshape_downsampled) %plot headshape
    ft_plot_sens(grad_trans, 'style', 'k*')
    ft_plot_vol(headmodel_singleshell,  'facecolor', 'cortex', 'edgecolor', 'none'); alpha 1; camlight
    view([90,0]); title('After Coreg');
    print([coreg_output 'headmodel_quality'],'-dpdf');

    fprintf('\nSaving the necessary data\n');
    save ([coreg_output 'headmodel_singleshell.mat'], 'headmodel_singleshell');
    save ([coreg_output 'mri_realigned_transformed.mat'], 'mri_realigned');
    save ([coreg_output 'trans_matrix.mat'], 'trans_matrix');
    save ([coreg_output 'grad_trans.mat'], 'grad_trans');
    save ([coreg_output 'matlab.mat']); % save all variables

    fprintf('\nCompleted MEMES.m - check the output for quality control\n');

    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Subfunctions

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [shape] = parsePolhemus(elpfile,hspfile)
        
        fid1 = fopen(elpfile);
        C = fscanf(fid1,'%c');
        fclose(fid1);
        
        E = regexprep(C,'\r','xx');
        E = regexprep(E,'\t','yy');
        
        returnsi = strfind(E,'xx');
        tabsi = strfind(E,'yy');
        sensornamesi = strfind(E,'%N');
        fiducialsstarti = strfind(E,'%F');
        lastfidendi = strfind(E(fiducialsstarti(3):fiducialsstarti(length(fiducialsstarti))+100),'xx');
        fiducialsendi = fiducialsstarti(1)+strfind(E(fiducialsstarti(1):fiducialsstarti(length(fiducialsstarti))+lastfidendi(1)),'xx');
        
        NASION = E(fiducialsstarti(1)+4:fiducialsendi(1)-2);
        NASION = regexprep(NASION,'yy','\t');
        NASION = str2num(NASION);
        
        LPA = E(fiducialsstarti(2)+4:fiducialsendi(2)-2);
        LPA = regexprep(LPA,'yy','\t');
        LPA = str2num(LPA);
        
        RPA = E(fiducialsstarti(3)+4:fiducialsendi(3)-2);
        RPA = regexprep(RPA,'yy','\t');
        RPA = str2num(RPA);
        
        LPAredstarti = strfind(E,'LPAred');
        LPAredendi = strfind(E(LPAredstarti(1):LPAredstarti(length(LPAredstarti))+45),'xx');
        LPAred = E(LPAredstarti(1)+11:LPAredstarti(1)+LPAredendi(2)-2);
        LPAred = regexprep(LPAred,'yy','\t');
        LPAred = str2num(LPAred);
        
        RPAyelstarti = strfind(E,'RPAyel');
        RPAyelendi = strfind(E(RPAyelstarti(1):RPAyelstarti(length(RPAyelstarti))+45),'xx');
        RPAyel = E(RPAyelstarti(1)+11:RPAyelstarti(1)+RPAyelendi(2)-2);
        RPAyel = regexprep(RPAyel,'yy','\t');
        RPAyel = str2num(RPAyel);
        
        PFbluestarti = strfind(E,'PFblue');
        PFblueendi = strfind(E(PFbluestarti(1):PFbluestarti(length(PFbluestarti))+45),'xx');
        PFblue = E(PFbluestarti(1)+11:PFbluestarti(1)+PFblueendi(2)-2);
        PFblue = regexprep(PFblue,'yy','\t');
        PFblue = str2num(PFblue);
        
        LPFwhstarti = strfind(E,'LPFwh');
        LPFwhendi = strfind(E(LPFwhstarti(1):LPFwhstarti(length(LPFwhstarti))+45),'xx');
        LPFwh = E(LPFwhstarti(1)+11:LPFwhstarti(1)+LPFwhendi(2)-2);
        LPFwh = regexprep(LPFwh,'yy','\t');
        LPFwh = str2num(LPFwh);
        
        RPFblackstarti = strfind(E,'RPFblack');
        RPFblackendi = strfind(E(RPFblackstarti(1):end),'xx');
        RPFblack = E(RPFblackstarti(1)+11:RPFblackstarti(1)+RPFblackendi(2)-2);
        RPFblack = regexprep(RPFblack,'yy','\t');
        RPFblack = str2num(RPFblack);
        
        allfids = [NASION;LPA;RPA;LPAred;RPAyel;PFblue;LPFwh;RPFblack];
        fidslabels = {'NASION';'LPA';'RPA';'LPAred';'RPAyel';'PFblue';'LPFwh';'RPFblack'};
        
        fid2 = fopen(hspfile);
        C = fscanf(fid2,'%c');
        fclose(fid2);
        E = regexprep(C,'\r','xx'); %replace returns with "xx"
        E = regexprep(E,'\t','yy'); %replace tabs with "yy"
        returnsi = strfind(E,'xx');
        tabsi = strfind(E,'yy');
        
        headshapestarti = strfind(E,'position of digitized points');
        headshapestartii = strfind(E(headshapestarti(1):end),'xx');
        headshape = E(headshapestarti(1)+headshapestartii(2)+2:end);
        headshape = regexprep(headshape,'yy','\t');
        headshape = regexprep(headshape,'xx','');
        headshape = str2num(headshape);
        
        shape.pnt = headshape;
        shape.fid.pnt = allfids;
        shape.fid.label = fidslabels;
        
        %convert to BESA style coordinates so can use the .pos file or sensor
        %config from .con
        shape.pnt = cat(2,fliplr(shape.pnt(:,1:2)),shape.pnt(:,3)).*1000;
        %shape.pnt = shape.pnt(1:length(shape.pnt)-15,:); % get rid of nose points may want to alter or comment this depending on your digitisation
        %shape.pnt = shape.pnt*1000;
        neg = shape.pnt(:,2)*-1;
        shape.pnt(:,2) = neg;
        
        shape.fid.pnt = cat(2,fliplr(shape.fid.pnt(:,1:2)),shape.fid.pnt(:,3)).*1000;
        %shape.fid.pnt = shape.fid.pnt*1000;
        neg2 = shape.fid.pnt(:,2)*-1;
        shape.fid.pnt(:,2) = neg2;
        shape.unit='mm';
        shape = ft_convert_units(shape,'cm');
        
        new_name2 = ['shape.mat'];
        save ([coreg_output new_name2], 'shape');
    end

    function [R,T,Yf,Err] = rot3dfit(X,Y)
        %ROT3DFIT Determine least-square rigid rotation and translation.
        % [R,T,Yf] = ROT3DFIT(X,Y) permforms a least-square fit for the
        % linear form
        %
        % Y = X*R + T
        %
        % where R is a 3 x 3 orthogonal rotation matrix, T is a 1 x 3
        % translation vector, and X and Y are 3D points sets defined as
        % N x 3 matrices. Yf is the best-fit matrix.
        %
        % See also SVD, NORM.
        %
        % rot3dfit: Frank Evans, NHLBI/NIH, 30 November 2001
        %
        
        % ROT3DFIT uses the method described by K. S. Arun, T. S. Huang,and
        % S. D. Blostein, "Least-Squares Fitting of Two 3-D Point Sets",
        % IEEE Transactions on Pattern Analysis and Machine Intelligence,
        % PAMI-9(5): 698 - 700, 1987.
        %
        % A better theoretical development is found in B. K. P. Horn,
        % H. M. Hilden, and S. Negahdaripour, "Closed-form solution of
        % absolute orientation using orthonormal matrices", Journal of the
        % Optical Society of America A, 5(7): 1127 - 1135, 1988.
        %
        % Special cases, e.g. colinear and coplanar points, are not
        % implemented.
        
        %error(nargchk(2,2,nargin));
        narginchk(2,2); %PFS Change to update
        if size(X,2) ~= 3, error('X must be N x 3'); end;
        if size(Y,2) ~= 3, error('Y must be N x 3'); end;
        if size(X,1) ~= size(Y,1), error('X and Y must be the same size'); end;
        
        % mean correct
        
        Xm = mean(X,1); X1 = X - ones(size(X,1),1)*Xm;
        Ym = mean(Y,1); Y1 = Y - ones(size(Y,1),1)*Ym;
        
        % calculate best rotation using algorithm 12.4.1 from
        % G. H. Golub and C. F. van Loan, "Matrix Computations"
        % 2nd Edition, Baltimore: Johns Hopkins, 1989, p. 582.
        
        XtY = (X1')*Y1;
        [U,S,V] = svd(XtY);
        R = U*(V');
        
        % solve for the translation vector
        
        T = Ym - Xm*R;
        
        % calculate fit points
        
        Yf = X*R + ones(size(X,1),1)*T;
        
        % calculate the error
        
        dY = Y - Yf;
        Err = norm(dY,'fro'); % must use Frobenius norm
    end

    function [output] = ft_transform_geometry_PFS_hacked(transform, input)
        
        % FT_TRANSFORM_GEOMETRY applies a homogeneous coordinate transformation to
        % a structure with geometric information, for example a volume conduction model
        % for the head, gradiometer of electrode structure containing EEG or MEG
        % sensor positions and MEG coil orientations, a head shape or a source model.
        %
        % The units in which the transformation matrix is expressed are assumed to
        % be the same units as the units in which the geometric object is
        % expressed. Depending on the input object, the homogeneous transformation
        % matrix should be limited to a rigid-body translation plus rotation
        % (MEG-gradiometer array), or to a rigid-body translation plus rotation
        % plus a global rescaling (volume conductor geometry).
        %
        % Use as
        %   output = ft_transform_geometry(transform, input)
        %
        % See also FT_WARP_APPLY, FT_HEADCOORDINATES
        
        % Copyright (C) 2011, Jan-Mathijs Schoffelen
        %
        % This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
        % for the documentation and details.
        %
        %    FieldTrip is free software: you can redistribute it and/or modify
        %    it under the terms of the GNU General Public License as published by
        %    the Free Software Foundation, either version 3 of the License, or
        %    (at your option) any later version.
        %
        %    FieldTrip is distributed in the hope that it will be useful,
        %    but WITHOUT ANY WARRANTY; without even the implied warranty of
        %    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
        %    GNU General Public License for more details.
        %
        %    You should have received a copy of the GNU General Public License
        %    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
        %
        % $Id: ft_transform_geometry.m$
        
        % flg rescaling check
        allowscaling = ~ft_senstype(input, 'meg');
        
        % determine the rotation matrix
        rotation = eye(4);
        rotation(1:3,1:3) = transform(1:3,1:3);
        
        if any(abs(transform(4,:)-[0 0 0 1])>100*eps)
            error('invalid transformation matrix');
        end
        
        %%### get rid of this accuracy checking below as some of the transformation
        %%matricies will be a bit hairy###
        if ~allowscaling
            % allow for some numerical imprecision
            %if abs(det(rotation)-1)>1e-6%100*eps
            %if abs(det(rotation)-1)>100*eps  % allow for some numerical imprecision
            %error('only a rigid body transformation without rescaling is allowed');
            %end
        end
        
        if allowscaling
            % FIXME build in a check for uniform rescaling probably do svd or so
            % FIXME insert check for nonuniform scaling, should give an error
        end
        
        tfields   = {'pos' 'pnt' 'o' 'coilpos' 'chanpos' 'chanposold' 'chanposorg' 'elecpos', 'nas', 'lpa', 'rpa', 'zpoint'}; % apply rotation plus translation
        rfields   = {'ori' 'nrm'     'coilori' 'chanori' 'chanoriold' 'chanoriorg'};                                          % only apply rotation
        mfields   = {'transform'};           % plain matrix multiplication
        recfields = {'fid' 'bnd' 'orig'};    % recurse into these fields
        % the field 'r' is not included here, because it applies to a volume
        % conductor model, and scaling is not allowed, so r will not change.
        
        fnames    = fieldnames(input);
        for k = 1:numel(fnames)
            if ~isempty(input.(fnames{k}))
                if any(strcmp(fnames{k}, tfields))
                    input.(fnames{k}) = apply(transform, input.(fnames{k}));
                elseif any(strcmp(fnames{k}, rfields))
                    input.(fnames{k}) = apply(rotation, input.(fnames{k}));
                elseif any(strcmp(fnames{k}, mfields))
                    input.(fnames{k}) = transform*input.(fnames{k});
                elseif any(strcmp(fnames{k}, recfields))
                    for j = 1:numel(input.(fnames{k}))
                        input.(fnames{k})(j) = ft_transform_geometry(transform, input.(fnames{k})(j));
                    end
                else
                    % do nothing
                end
            end
        end
        output = input;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SUBFUNCTION that applies the homogeneous transformation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [new] = apply(transform, old)
        old(:,4) = 1;
        new = old * transform';
        new = new(:,1:3);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % rotate_about_z - make a rotation matix for arbitrary rotation in degrees
    % around z axis
    %
    % Written by Paul Sowman Oct 2017 (http://web.iitd.ac.in/~hegde/cad/lecture/L6_3dtrans.pdf - page 4)
    %
    % INPUTS:
    % - deg        = degrees of rotation required
    %
    % OUTPUTS:
    % - rmatx      = a 4*4 rotation matrix for deg degrees about z
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function rmatx = rotate_about_z(deg)
        
        rad   = pi/180 * deg; %deg2rad(deg); 
        rmatx = [cos(rad) sin(rad) 0 0;-sin(rad) cos(rad) 0 0;0 0 1 0;0 0 0 1];
    end

    function [headshape_downsampled] = downsample_headshape_noface(path_to_headshape,numvertices,sensors)
        % Get headshape
        headshape = ft_read_headshape(path_to_headshape);
        % Convert to cm
        headshape = ft_convert_units(headshape,'cm');
        % Convert to BESA co-ordinates
        headshape.pos = cat(2,fliplr(headshape.pos(:,1:2)),headshape.pos(:,3));
        headshape.pos(:,2) = headshape.pos(:,2).*-1;
        
        % Get indices of facial points (up to 4cm above nasion)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Is 4cm the correct distance?
        % Possibly different for child system?
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        count_facialpoints = find(headshape.pos(:,3)<4);
        if isempty(count_facialpoints)
            disp('CANNOT FIND ANY FACIAL POINTS');
        else
            facialpoints = headshape.pos(count_facialpoints,:,:);
            rrr = 1:4:length(facialpoints);
            facialpoints = facialpoints(rrr,:); clear rrr;
        end
        
        % Remove facial points for now
        headshape.pos(count_facialpoints,:) = [];
        
        % Create mesh out of headshape downsampled to x points specified in the
        % function call
        cfg.numvertices = numvertices;
        cfg.method = 'headshape';
        cfg.headshape = headshape.pos;
        mesh = ft_prepare_mesh(cfg, headshape);
        
        % Replace the headshape info with the mesh points
        headshape.pos = mesh.pos;
        
        % Create figure for quality checking
        figure; subplot(2,2,1);ft_plot_mesh(mesh); hold on;
        title('Downsampled Mesh');
        view(0,0);
        subplot(2,2,2);ft_plot_mesh(headshape); hold on;
        title('Downsampled Headshape View 1');
        view(0,0);
        subplot(2,2,3);ft_plot_mesh(headshape); hold on;
        title('Downsampled Headshape View 2');
        view(90,0);
        subplot(2,2,4);ft_plot_mesh(headshape); hold on;
        title('Downsampled Headshape View 3');
        view(180,0);
        print([coreg_output 'headshape_quality'],'-dpdf');
        
        % Add in names of the fiducials from the sensor
        headshape.fid.label = {'NASION','LPA','RPA'};
        
        % Convert fiducial points to BESA
        headshape.fid.pos = cat(2,fliplr(headshape.fid.pos(:,1:2)),headshape.fid.pos(:,3));
        headshape.fid.pos(:,2) = headshape.fid.pos(:,2).*-1;
        
        % Plot for quality checking
        figure;ft_plot_sens(sensors) %plot channel position : between the 1st and 2nd coils
        ft_plot_headshape(headshape) %plot headshape
        view(0,0);
        print([coreg_output 'headshape_quality2'],'-dpdf');
        
        % Export filename
        headshape_downsampled = headshape;
        
    end

end




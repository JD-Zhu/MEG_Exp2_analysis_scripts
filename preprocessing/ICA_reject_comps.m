% Manually select which ICA components to reject:
% eye blinks, other large muscle artefacts (e.g. jaw clenching, hand movements)

function [data_clean] = ICA_reject_comps(data, comp, lay)

    % Plot the ICA components, so we can identify which comps to remove

    % change the colourmap
    ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
    colormap(flipud(brewermap(64, 'RdBu'))) % change the colormap
    
    cfg = [];
    cfg.component = 1:49; % this fits on 1 page      
    cfg.layout = lay;
    cfg.marker = 'off';
    cfg.comment = 'no';
    ft_topoplotIC (cfg, comp);

    cfg          = [];
    %cfg.channel  = 1:10; % components to be plotted
    cfg.viewmode = 'component';
    cfg.blocksize = 60; % display 60-sec segments
    cfg.layout   = lay;
    ft_databrowser(cfg, comp);
    drawnow; pause;
    
    
    % Remove components
    diary on;
    success = false;
    
    while ~success
        % Ask user to specify the components to be removed
        disp('Enter components in the form [1 2 3]'); drawnow;
        comps2remove = input('Which components would you like to remove?\n'); drawnow;

        % Remove these components
        cfg           = [];
        cfg.component = comps2remove; %these are the components to be removed
        data_clean    = ft_rejectcomponent(cfg, comp, data);


        % Quality check: plot the continuous data BEFORE and AFTER rejecting comps
        cfg          = [];
        cfg.channel  = 'all'; % components to be plotted
        cfg.viewmode = 'vertical';
        cfg.blocksize = 60; % display 60-sec segments
        cfg.ylim     = [ -5e-13 5e-13 ];
        cfg.layout   = lay;

        % plot BEFORE
        ft_databrowser(cfg, data);
        % plot AFTER
        ft_databrowser(cfg, data_clean);
        drawnow;

        
        % Let user decide whether they want to retry with diff comps
        prompt = ['\nCompare the datasets with and without IC removal:\n' ...
            'If it is OK, press Y to continue.\nIf not ok, press N, ' ...
            'then you will be asked to select the components again.\n\n'];

        answer = input(prompt, 's'); % read keyboard input as a string
        if strcmp(answer, 'N')
            success = false;
        else % for any other key press, we treat it as 'Yes'
            success = true;
        end
    end
    
    diary off;

end
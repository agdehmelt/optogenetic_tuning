list = {'Figure5B','Figure5C','Figure5E'};
[indx,tf] = listdlg('PromptString','Generate simulations related to:','SelectionMode','single','ListString',list);

figure=char(list(indx));

switch figure
    case {'Figure5B'}
        settings='single_low_gt';
        ca_sim_RGM;
    case {'Figure5C'}
        settings='scan_gefpm';
        ca_sim_RGM;
        settings='scan_gefcyto';
        ca_sim_RGM;
        settings='scan_myopm';
        ca_sim_RGM;
        settings='scan_myocyto';
        ca_sim_RGM;
        settings='scan_rhopm';
        ca_sim_RGM;
        settings='scan_rhocyto';
        ca_sim_RGM;
    case {'Figure5E'}
        settings='single_high_gt';
        ca_sim_RGM;
end

        
    
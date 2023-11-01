function Spock_CombineStacksBVcorrectTrial(folder_path,save_fn,parameter_class)
% folder_path = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA008/DA008_101823/DA008_101823_img/DA008_101823_1'; 
% save_fn = '/Volumes/buschman/Rodent Data/Behavioral_dynamics_cj/DA008/DA008_101823/DA008_101823_img/DA008_101823_1/DA008_101823_1_dff_combined.mat';
% parameter_class = 'general_params_example';  

    if ~ispc
        addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/Widefield_Imaging_Analysis/'));
    end
    try
        [stack, opts] = CombineStacks(folder_path,parameter_class);       
        %perform hemodynamic correction and make dff
        if numel(unique(opts.wavelength_pattern))>1 %if multiple wavelengths used
           [dff, dff_b, ~] = HemodynamicCorrectionBVcorrectTrial(stack, opts); 
           %ImpactOfHemoCorrection(dff,dff_b,dff_h)
        else
           fprintf('\n No hemodynamic correction');
           dff_b = makeDFF(stack, opts); 
           dff = [];
        end
        
        fprintf('\n Saving data');
        %Save off corrected data if available
        if ~isempty(dff)
            %save dff
            save(save_fn,'dff','opts','-v7.3');
        end
            
        %save off the uncorrected if desired
        if opts.save_uncorrected
            [path, fn] = fileparts(save_fn);
            dff = dff_b; %need to have it still named dff. 
            fn = [path filesep fn '_dff_uncorrected.mat'];
            save(fn,'dff','opts','-v7.3');
        end
        
%         %make figures
%         fh = DetectDictalEvents(save_fn);
    catch
        gp = general_params; 
        warning('\nNo dffs with suffix %s were combined in %s',gp.dff_suffix, folder_path);
    end
end


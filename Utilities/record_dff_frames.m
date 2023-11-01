function record_dff_frames(tbytDat, filePathTrials, varargin)
% e.g., dff = tbytDat(21).dffsm
%By default this function iterate through all trials and record frames of
% each trial as png files. It can also run specific trials designated as varargin.  
if isempty(varargin)
    trials = 1:length(tbytDat); 
else
    trials = [varargin{1}(:)]'; 
end

%% record frames first
close all;
for t = trials % 1:length(tbytDat) % trial
    if ~isempty(tbytDat(t).dffsm)
        pngSubDir = fullfile(filePathTrials, sprintf('block_%d_trial_%d', tbytDat(t).cmosExpTrainI, t));

        if exist(pngSubDir, 'dir') == 7

        else
            mkdir(pngSubDir);
            for i = 1:size(tbytDat(t).dffsm,3)
                pngName = sprintf('frame_%d.png', i);
                figHandle = imageFrameWithNaNs(tbytDat(t).dffsm(:, :, i), [-2 2]); hold on;

                % Turn off the axes
                axis off;
                % Set the figure's background to white
                set(gcf, 'Color', 'w');

                if isfield(tbytDat, 'frameStimI')
                    if ~isempty(tbytDat(t).frameStimI) && tbytDat(t).frameStimI(i) % for visual trials
                        if tbytDat(t).evtType == 1 % common visual stim (draw 45 degree lines at the upper left corner)
                            insertgrating45(figHandle, tbytDat(t).dffsm(:, :, i))
                        elseif tbytDat(t).evtType == 2 % common visual stim (draw 135 degree lines at the upper left corner)
                            insertgrating135(figHandle, tbytDat(t).dffsm(:, :, i))
                        end
                        hold off;
                    end
                end
                % Capture the current figure with dots
                frameLabled = getframe(gca);
                frameLabled = frameLabled.cdata;

                % save the figure with dots

                imwrite(frameLabled, fullfile(pngSubDir, pngName));

                fprintf('Frame #%d of trial #%d is labeled and saved.\n', i, t);
                close all
            end
        end
    end
end

end

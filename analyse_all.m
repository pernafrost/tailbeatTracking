
clear all;
close all;

addpath('/Applications/dip/common/dipimage')
dip_initialise

pathToSeries = [pwd, filesep, 'input_images/test_series/'];

allDirs = dir(pathToSeries);

for dd = 1 : length(allDirs)
    currentDir = [pathToSeries, allDirs(dd).name];
    if exist(currentDir, 'dir') && currentDir(end) ~= '.' && currentDir(1) ~= '.'

        
        file = fullfile(currentDir, '00000001.png')
        [pathName, trialName, extensionName] = fileparts(currentDir);
        trialName = [trialName, extensionName];
        disp(sprintf('Processing now "%s"', trialName));
        medianImageInput = 'none';
        medianImageOutput = 'none';
                
        % I don't use the median image because the background is
        % saturated!!!
%         medianImageInput=[pathName, filesep, 'median_image_', trialName(8:end), '.mat'];
%         medianImageOutput = medianImageInput
%         if ~exist(medianImageInput, 'file')
%             medianImageInput = 'none'
%         end
        
        resultsFile = [pathName, filesep, 'results_', trialName(8:end), '.mat']
        body_synchronization_analysis(file, medianImageInput, medianImageOutput, resultsFile);
        
        check_tracking_results(resultsFile, file, strrep(currentDir, 'images_', 'results_'));
    end
end



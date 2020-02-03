function [medianImage]=compute_median_image(inputFilePath, fileList, useChannel, nStepsForMedianImage, varargin)

if size(varargin, 2) > 0
    fileNameForMedianImageIn = varargin{1};
else
    fileNameForMedianImageIn = '';
end

if size(varargin, 2) > 1
    fileNameForMedianImageOut = varargin{2};
else
    fileNameForMedianImageOut = '';
end

if ~isempty(fileNameForMedianImageIn) % if I already input a file name from which to read the median image
%     [medianImagePath, medianImageFile, medianImageExtension,~] = fileparts(fileNameForMedianImageIn);
%    	medianImagePath = [medianImagePath, filesep];
%     medianImageFile = [medianImageFile, medianImageExtension];
%     medianImageFile
%     medianImagePath
%     fullfile(medianImagePath, medianImageFile)
%     exist(medianImageFile,'file')
%     exist(medianImagePath,'dir')
    
    if ~exist(fileNameForMedianImageIn,'file')
        disp('no file selected; switching to computing median image from the image sequence')
        medianImage= compute_median_image_from_images(inputFilePath, fileList, useChannel, nStepsForMedianImage, fileNameForMedianImageOut);
    else
        load(fileNameForMedianImageIn);
    end
    
else % if I do not input a file name, ask the user
    medianImageDlg = questdlg('Do you want to read the median image from file?', ...
        'Median image question', ...
        'Yes', 'No', 'No');
    switch medianImageDlg,
        case 'Yes',
            [medianImageFile, medianImagePath] = uigetfile('*.mat;*.MAT', 'Pick the file with the median image');
            medianImageFile = fullfile(medianImagePath, medianImageFile);
            if ~exist(medianImageFile,'file')
                disp('no file selected; switching to computing median image from the image sequence')
                medianImage= compute_median_image_from_images(inputFilePath, fileList, useChannel, nStepsForMedianImage, fileNameForMedianImageOut);
            else
                load(medianImageFile);
            end
        case 'No',
            medianImage = compute_median_image_from_images(inputFilePath, fileList, useChannel, nStepsForMedianImage, fileNameForMedianImageOut);
    end % switch
    
end;






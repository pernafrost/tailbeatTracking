%% Program body synchronization analysis

% the argument can be the complete file path and name
function [] = body_synchronization_analysis(varargin)
% clear all;
close all;

if size(varargin, 2) > 0
    completeFileNameAndPath = varargin{1};
else
    completeFileNameAndPath = '';
end

if size(varargin, 2) > 1
    fileNameForMedianImageIn = varargin{2};
else
    fileNameForMedianImageIn = 'none';
end

if size(varargin, 2) > 2
    fileNameForMedianImageOut = varargin{3};
else
    fileNameForMedianImageOut = '';
end

if size(varargin, 2) > 3
    fileNameForResults = varargin{4};
else
    fileNameForResults = '';
end
    
%    completeFileNameAndPath = []; fileNameForMedianImageIn = '';
%   fileNameForMedianImageOut = ''; fileNameForResults = '';

showFigures = 0; % 0=no, 1= only big figure, 2=all figures
useChannel = 1; % 0=red, 1=green, 2=blue % channel was 1 in the image with the large basin


%% Input dir with the images
% File names should be sorted according to the time
% I ask the user to pick a file of the sequence and then I will process all
% the files with the same name (except for the final numbering).

if isempty(completeFileNameAndPath)
    [inputFileName, inputFilePath]=uigetfile('*.PNG;*.png;*.jpg;*.JPG', 'Pick one image of the sequence');
else
    [inputFilePath, inputFileName, inputFileExtension] = fileparts(completeFileNameAndPath);
    inputFilePath = [inputFilePath, filesep];
    inputFileName = [inputFileName, inputFileExtension];
end


if ~exist(fullfile(inputFilePath, inputFileName),'file')
    disp('No valid file name');
    return;
end

imgNumPosition = regexp(inputFileName(1:end-4), '\D');
if isempty(imgNumPosition) % there are only numbers
    baseFileName ='';
else
    baseFileName = inputFileName(1:imgNumPosition(end));
end
[~,~,inputFileExtension] = fileparts(inputFileName);


fileList = dir([inputFilePath, baseFileName,'*', inputFileExtension]);



% initialize dip library

start_dip();

% compute median image
if strcmp(fileNameForMedianImageIn, 'none') % if the name of the median image file is none, we do not compute the median image
    medianImage = [];
else
    nStepsForMedianImage = 40;
    medianImage = compute_median_image(inputFilePath, fileList, useChannel, nStepsForMedianImage, fileNameForMedianImageIn, fileNameForMedianImageOut);
end


results = cell(length(fileList), 1);

labelCounter = 0;

for jj = 1:length(fileList)
    fileNameAndPath = fullfile(inputFilePath, fileList(jj).name);
    disp(sprintf('processing now %s', fileList(jj).name));
    img = imread(fileNameAndPath);
    imgDip = dip_image(img);
    % try

    [msr, componentLines, rotAngle, bodyCoords, fitParams, objectLabels, labelCounter] = process_image(imgDip, useChannel, medianImage, labelCounter, showFigures);
    results{jj}.msr = msr;
    results{jj}.componentLines = componentLines;
    results{jj}.rotAngle = rotAngle;
    results{jj}.bodyCoords = bodyCoords;
    results{jj}.fitParams = fitParams;
    results{jj}.objectLabels = objectLabels;
    results{jj}.inputImage = fileNameAndPath;
%     objectLabels
%     labelCounter
    % end
    
%     feedBack = input('Is everything ok? ', 's');
%     if (feedBack == 'y')
%         close all;
%     else
%         disp('error...');
%         close all;
%     end
end

if isempty(fileNameForResults)
    [resultsFileName, resultsFilePath]=uiputfile('*.MAT;*.mat', 'Save results as');
    if isequal(resultsFileName,0) || isequal(resultsFilePath,0)
    disp('No valid file name');
    return;
    else
    fileNameForResults = fullfile(resultsFilePath, resultsFileName);
    end
end

save(fileNameForResults, 'results');






function [] = check_tracking_results(varargin)

size(varargin,2)


if size(varargin, 2) > 0 && exist(varargin{1}, 'file')
    resultsFile = varargin{1}
else
    [resultsFile, resultsDir] = uigetfile('*.mat', 'Pick the file with results');
    if isequal(resultsFile,0)
        return;
    end
    resultsFile = fullfile(resultsDir, filesep, resultsFile)
end


if size(varargin, 2) > 1 && exist(varargin{2}, 'file')
    imageFile = varargin{2};
    [inputFilePath, inputFileName, inputFileExtension] = fileparts(imageFile);
    inputFilePath = [inputFilePath, filesep]
    inputFileName = [inputFileName, inputFileExtension]
else
    [inputFileName, inputFilePath]=uigetfile('*.PNG;*.png;*.jpg;*.JPG', 'Pick one image of the sequence');
    if isequal(inputFileName,0) || isequal(inputFilePath,0)
        return;
    end
end

if size(varargin, 2) > 2 && exist(varargin{3}, 'dir')
    outputFolder = varargin{3};
else
    outputFolder = strrep(inputFilePath, 'input_images', 'result_images');
end


load(resultsFile, 'results');

isSavingTrjectoriesAsTxt = 1;
if isSavingTrjectoriesAsTxt
    fileTrajectories = fopen('trajectories.txt', 'w');
end

imgNumPosition = regexp(inputFileName(1:end-4), '\D');
if isempty(imgNumPosition) % there are only numbers
    baseFileName ='';
else
    baseFileName = inputFileName(1:imgNumPosition(end));
end
[~,~,inputFileExtension] = fileparts(inputFileName);


fileList = dir([inputFilePath, baseFileName,'*', inputFileExtension]);

start_dip();
myColours = 'rgbmycwk';
myLabelColours = copper(256);

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

for jj = 1:length(fileList)
    fileNameAndPath = fullfile(inputFilePath, fileList(jj).name);
    disp(sprintf('processing now %s', fileList(jj).name));
    img = imread(fileNameAndPath);
    imgDip = dip_image(img);
    imDisplay = dipshow(imgDip);
    diptruesize(imDisplay, 100);
    hold on;
    for ll = 1:length(results{jj}.componentLines)
         plot(results{jj}.componentLines{ll}(:,1),results{jj}.componentLines{ll}(:,2),myColours(mod(results{jj}.objectLabels(ll)-1, length(myColours))+1));
         
         
     
%        text(results{jj}.msr.data{1}(1,ll)-10, results{jj}.msr.data{1}(2,ll)+15,num2str(results{jj}.objectLabels(ll)) ,'color', myColours(mod(results{jj}.objectLabels(ll)-1, length(myColours))+1));
%        plot(results{jj}.msr.data{1}(1,ll), results{jj}.msr.data{1}(2,ll), 'Marker', 'o', 'MarkerSize', 2, 'Color', myColours(mod(results{jj}.objectLabels(ll)-1, length(myColours))+1)); 
         text(results{jj}.msr.Center(1,ll)-10, results{jj}.msr.Center(2,ll)+15,num2str(results{jj}.objectLabels(ll)) ,'color', myColours(mod(results{jj}.objectLabels(ll)-1, length(myColours))+1));
         plot(results{jj}.msr.Center(1,ll), results{jj}.msr.Center(2,ll), 'Marker', 'o', 'MarkerSize', 2, 'Color', myColours(mod(results{jj}.objectLabels(ll)-1, length(myColours))+1));
         
        if isSavingTrjectoriesAsTxt
         fprintf(fileTrajectories, '%d %03f %03f %d\n', jj, results{jj}.msr.Center(1,ll), results{jj}.msr.Center(2,ll), results{jj}.objectLabels(ll));
        end
        
    end
    
   print(gcf, '-dpng', fullfile(outputFolder,fileList(jj).name))
    
    %pause;
    close all;
end


if isSavingTrjectoriesAsTxt
fclose(fileTrajectories)
end



function medianImage = compute_median_image_from_images(inputFilePath, fileList, useChannel, nStepsForMedianImage, varargin)


if size(varargin, 2) > 0
    fileNameForMedianImageOut = varargin{1};
else
    fileNameForMedianImageOut = '';
end

clear imageArray;


disp('computing the median image may take quite a lot of time');
counter = 1;
if nStepsForMedianImage > length(fileList)
    allTheValues = 1:length(fileList);
else
    allTheValues = round(linspace(1, length(fileList), nStepsForMedianImage));
end

for jj = allTheValues
    % jj
    fileNameAndPath = fullfile(inputFilePath, fileList(jj).name);
    fprintf('.');
    img = imread(fileNameAndPath);
    imgGray = squeeze(img(:,:,useChannel));
    if (counter == 1)
        imageArray = zeros([size(imgGray), length(linspace(1, length(fileList), nStepsForMedianImage))]);
    end
    imageArray(:,:, counter) = imgGray;
    counter = counter + 1;
end

fprintf('\nfinished reading the images, now computing the median\nbe patient\n');

if size(imageArray,3) == 1
    medianImage = squeeze(imageArray(:,:,1));
else
    
    medianImage = squeeze(median(dip_image(imageArray),[],3));
end



if ~isempty(fileNameForMedianImageOut) % if I already input a file name from which to read the median image
    if strmatch(fileNameForMedianImageOut, 'none', 'exact')
        return;
    else
        save(fileNameForMedianImageOut, 'medianImage');
    end
    
else % if the file name was not passed in input, ask the user.
    
    medianImageDlg = questdlg('Do you want to save the median image to file?', ...
        'Median image question', ...
        'Yes', 'No', 'No');
    
    switch medianImageDlg,
        case 'Yes',
            [medianImageFile, medianImagePath] = uiputfile('*.mat;*.MAT', 'Select a file for the median image', fullfile(inputFilePath, 'median_image.mat'));
            if isequal(medianImageFile,0) || isequal(medianImagePath,0)
                disp('no file selected');
                disp('won''t save the median image');
                return;
            else
                save(medianImageFile, 'medianImage');
            end
        case 'No',
            disp('won''t save the median image')
            return;
    end % switch
    
    
end


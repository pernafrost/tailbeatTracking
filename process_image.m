function [msr, componentLines, rotAngle, bodyCoords, fitParams, objectLabels, labelCounter] = process_image(imgDip, useChannel, medianImage, labelCounter, varargin)
% disp('processing now');

persistent previousObjectLabels;
persistent previousObjectCoordinates;
% persistent labelCounter; % I prefer to pass labelCounter to this
% function, so that if I need to stop and restart the program I keep more
% control on the counter
% if isempty(labelCounter)
%     labelCounter = 0;
% end

thresholdDistance = 300;
optargin = size(varargin,2);
if optargin >= 1
    showFigures = varargin{1};
else
    showFigures = 1;
end

% imgDip

% red = squeeze(imgDip(:,:,0));
% green = squeeze(imgDip(:,:,1));
% blue = squeeze(imgDip(:,:,2));
%
% diphist(red);
% diphist(green);
% diphist(blue);

% binaryImg = threshold(gaussf(red, 3, 'best'), 'background', 15)
% binaryImg = threshold(gaussf(green, 3, 'best'),'triangle')
% binaryImg = threshold(gaussf(blue, 3, 'best'),'triangle')

% if I use a median image 

myGrayImage = squeeze(imgDip(:,:,useChannel));

if (sum(size(medianImage) == size(myGrayImage)) == 2) % if the medianImage and the gray image have the same size
    myGrayImage = medianImage - myGrayImage; % I use + sign because myGrayImage is already negative (I want the fish into the components)
end

% % test different thresholding methods
% threshold(gaussf(myGrayImage,1,'best'), 'minerror')
% threshold(gaussf(myGrayImage,1,'best'), 'background')
% threshold(gaussf(myGrayImage,1,'best'), 'isodata')
% threshold(gaussf(myGrayImage,1,'best'), 'otsu')
% threshold(gaussf(myGrayImage,1,'best'), 'triangle')
% threshold(gaussf(myGrayImage,1,'best'), 'background',3)
% pause;
% close all;

mySmoothGrayImage = gaussf(myGrayImage,1,'best');

if (sum(size(medianImage) == size(myGrayImage)) == 2) % if the medianImage and the gray image have the same size
    binaryImg = threshold(mySmoothGrayImage, 'fixed',40); % it was 40 for the video in the large basin (on the green channel), which was a bit darker
else
    
    mySmoothGrayImageClosed = closing(mySmoothGrayImage,50, 'rectangular') - mySmoothGrayImage; % rectangular filter is faster though a little bit less accurate
    mySmoothGrayImageReconstructed = reconstruction(mySmoothGrayImageClosed*(mySmoothGrayImageClosed > 100), mySmoothGrayImageClosed);
    binaryImg = mySmoothGrayImageReconstructed > 15; % threshold(mySmoothGrayImageReconstructed, 'fixed', 15);
    myGrayImage = mySmoothGrayImageReconstructed;
end

% binaryImg = threshold(mySmoothGrayImage, 'triangle')
labeledImage = label(binaryImg, inf, 500);
% labeledImage ~= 0
% diptruesize(gcf, 50);

% msr = measure(labeledImage,myGrayImage,{'GreyMu'});
try % I may get a problem with the complete measurement
    msr = measure(labeledImage,myGrayImage,{'Size', 'Perimeter', 'Minimum', 'Maximum', 'GreyMu', 'Gravity', 'Center', 'P2A', 'MajorAxes', 'PodczeckShapes', 'DimensionsEllipsoid', 'Feret', 'GreyMajorAxes', 'StdDev', 'Skewness', 'ExcessKurtosis'});
catch
    % msr = [];
    % msr = measure(labeledImage,myGrayImage,{'Size', 'Perimeter', 'Minimum', 'Maximum', 'GreyMu', 'Gravity', 'Center', 'P2A', 'MajorAxes', 'PodczeckShapes', 'DimensionsEllipsoid', 'Feret'});
    msr = measure(labeledImage,myGrayImage,{'Size', 'Minimum', 'Maximum', 'GreyMu', 'Gravity', 'Center', 'P2A'});
end
imsz = imsize(binaryImg);

componentLines=cell(length(msr), 1);
bodyCoords=cell(length(msr), 1);
rotAngle=zeros(length(msr),1);
fitParams=cell(length(msr),1);
objectLabels=zeros(length(msr),1);

fishBodyLength=zeros(length(msr),1);

% exclude coordinates close to the border?
if ~isempty(msr)
    
    
    
    
    for ll = 1:length(msr)
        
        %        disp(['Fish number ',num2str(ll)])
        %        disp(['   position: (',num2str(msr(ll).Gravity(1),'%.3f'),',', num2str(msr(ll).Gravity(2),'%.3f'),')']);
        
        %         if ((msr(ll).minimum(1) == 0) || (msr(ll).minimum(2) == 0) || (msr(ll).maximum(1) == (imsz(1) -1)) || (msr(ll).maximum(2) == (imsz(2) -1)))
        %             continue; % the fish is cut by the border
        %         end
        
        
        % Cut out a portion of the image around the fish.
        cc = [max(msr(ll).Minimum'-2,0); min(msr(ll).Maximum'+2,imsz -1)];
        sz = diff(cc)+1;
        if mod(sz(1),2)==0 % Even in size along x, extend by 1 pixel.
            if cc(1,1)==0 && cc(2,1) == imsz(1) -1 % if it is as large as the whole image reduce the size
                cc(2,1) = cc(2,1) -1;
            elseif cc(1,1)==0 % if it starts at zero, increase size at end
                cc(2,1) = cc(2,1)+1;
            else % if it does not start at zero, increase size at lower point
                cc(1,1) = cc(1,1)-1;
            end
        end
        
        if mod(sz(2),2)==0 % Even in size along y, extend by 1 pixel.
            if cc(1,2)==0 && cc(2,2) == imsz(2) -1
                cc(2,2) = cc(2,2) -1;
            elseif cc(1,2) == 0
                cc(2,2) = cc(2,2)+1;
            else
                cc(1,2) = cc(1,2)-1;
            end
        end
        
        
        
        bin = labeledImage(cc(1,1):cc(2,1),cc(1,2):cc(2,2)) == ll;
        gv  = myGrayImage(cc(1,1):cc(2,1),cc(1,2):cc(2,2));
        
        % Note: we want odd patch size just so that the rotation
        % later is easier to understand and undo.
        
        % Determine the main orientation of the fish.
        t = msr(ll).GreyMu;
        % The measure GreyMu uses the greyvalues within the object's
        % mask, so for this type of object, where the background
        % intensity is +/- 0, it is very much independent of a good
        % segmentation. It returns the values of the inertia tensor,
        % from which we compute the orientation unafected by
        % discretization.
        
        t = [t(1),t(2);t(2),t(3)];
        [v,e] = eig(t);
        if e(1,1)<e(2,2)
            v = v(1,:);
        else
            v = v(2,:);
        end
        rotAngle(ll) = - atan2(v(2),v(1)); % rotation angle phi
        %       disp(['   orientation: ',num2str(rotAngle(ll)),' rad'])
        % This won't always put the head at the same location. We take
        % care of that later.
        
        % Rotate the image
        gvRot = rotation(gv,rotAngle(ll),3,'3-cubic','zero');
        binRot = rotation(+bin,rotAngle(ll),3,'zoh','zero');
        
        gvRotSmooth = gaussf(gvRot, 5, 'best');
        
        % Make a new mask to find out which lines to measure on.
        % binRot2 = threshold(gvRotSmooth) | binRot; % QUESTION: should I use | or & ? Or I can just use binRot?
        binRot2 = binRot ~= 0;
        x = xx(binRot2,'corner');
        x = x(binRot2);
        range = [min(x),max(x)]; % this is the min and max coords in object
        
        
        
        % Determine position along each of the lines.
        xCoords = range(1):range(2);
        yCoords = zeros(size(xCoords));
        value = yCoords;
        for ii=1:length(xCoords)
            % For each integer x coordinate, we extract a 1D image and determine
            % the location of the maximum with sub-pixel accuracy.
            line = squeeze(gvRotSmooth(xCoords(ii),:));
            mask = squeeze(binRot2(xCoords(ii),:));
            [value(ii),position] = max(line,mask);
            yCoords(ii) = subpixlocation(line,position,'parabolic','maximum');
        end
        
        
        
        % Determine which side the head is on.
        n = floor(length(value)/2);
        p = mean(value(1:n)) - mean(value(n:end));
        % Even better would be to look at smaller patches, maybe
        % 1/5 to 1/3 and 2/3 to 4/5 of the length.
        % Maybe this:
        % p = p(3:end-2);
        % n = ceil(length(p)/3);
        % p = mean(p(1:n)) - mean(p(2*n:end));
        invertedRotAngle = 0;
        if p > 10
            % Head on the left: don't do anything.
            %           disp(['   direction: ',num2str(rotAngle(ll)),' rad'])
            invertedRotAngle = 0;
        elseif p < -10
            % Head on the right: flip.
            xCoords = xCoords(end:-1:1);
            yCoords = yCoords(end:-1:1);
            %            disp(['   direction: ',num2str(-rotAngle(ll)),' rad'])
            invertedRotAngle = 1;
        else
            % Difference between head and tail is small.
            %            disp('   couldn''t determine direction!')
            invertedRotAngle = 0;
            % I don't know if the +/- 5 grey values is a good threshold
            % to determine uncertainty or not. We'd need to look at lots
            % of fish to be sure about noise levels and so on.
        end
        
        
        bodyCoords{ll}.x = xCoords;
        bodyCoords{ll}.y = yCoords;
        
        
        %       fishBodyLength(ll) = curve_length(xCoords, yCoords);
        
        
        
        % Fit some function on the fish body
        % options.MaxFunEvals = 20000;
        n1 = floor(length(value)*3/4); % I use the value n1 for selecting the portion over which to fit the curve
        data = [xCoords(n1+1:end)', yCoords(n1+1:end)'];
        
        %         Q = @(params) fitting_my_function(params, data);
        %         params_init = [(rand(1)*2 -1), exp(mvnrnd([log(5), log(0.9), log(1.5), log(1)], eye(4)))]
        %         params_best = fminsearch(Q, params_init)
        %
        %         fitParams{ll}.phi = params_best(1);
        %         fitParams{ll}.omega = params_best(2);
        %         fitParams{ll}.A = params_best(3);
        %         fitParams{ll}.gamma = params_best(4);
        %         fitParams{ll}.beta = params_best(5);
        baseLine = nanmean(yCoords);
        %        fitParams{ll}.zeroCrossing = find(diff(sign(data(:,2)-baseLine)), 1, 'last');
        %        fitParams{ll}.meanDeviationFromTailToZeroCrossing = mean(data(fitParams{ll}.zeroCrossing+1:end, 2) -baseLine);
        %        [fitParams{ll}.maximum, fitParams{ll}.posMaximum] = max(data(:,2)-baseLine);
        %        [fitParams{ll}.minimum, fitParams{ll}.posMinimum] = min(data(:,2)-baseLine);
        
        fitParams{ll}.meanDeviation = mean(data(:, 2) - baseLine);
        fitParams{ll}.tailFinPosition = data(end, 2) - baseLine;
        segmentLength = 10; % a segment is 10 pixels
        for jj = 1:7
            try
                fitParams{ll}.lastSegment(1:2, jj) = polyfit(yCoords(end - segmentLength*(jj)+1:end - segmentLength*(jj-1), 1), yCoords(end-segmentLength*(jj)+1:end - segmentLength*(jj-1), 2),1);
            catch
                fitParams{ll}.lastSegment(1:2,jj) = [NaN, NaN];
            end
        end
        
        
        % Display, to see if it works.
        if showFigures >= 2
            hSmallImg = dipshow(gvRot);
            diptruesize(hSmallImg, 500);
            hold on
            plot(xCoords,yCoords,'r');
            plot(xCoords,mean(yCoords)*ones(size(xCoords)), 'w');
            plot(xCoords(1), yCoords(1),'or');
            plot(data(1,1), data(1, 2),'ow');
            % plot(data(end-19:end-10,1), 135.8363 + 0.3563*data(end-19:end-10,1), 'b');
            %            plot(data(fitParams{ll}.zeroCrossing,1), data(fitParams{ll}.zeroCrossing, 2),'ow');
            %            plot(data(fitParams{ll}.zeroCrossing+1:end, 1), data(fitParams{ll}.zeroCrossing+1:end, 2), 'g');
            %            t = pi*(data(:,1) - min(data(:,1))) / (max(data(:,1)) - min(data(:,1)));
            %            plot(data(:,1), mean(data(:,2)) + (max(data(:,1)) - min(data(:,1))) * (fitParams{ll}.A * cos(fitParams{ll}.omega*t.^fitParams{ll}.gamma + fitParams{ll}.phi)) .* exp(-fitParams{ll}.beta*t), 'g');
        end
        
        
        
        % Bring back found coordinates to image coordinates.
        fitParams{ll}.orig = floor(imsize(gvRotSmooth)/2);
        fitParams{ll}.cc = cc;
        componentLines{ll} = [xCoords' - fitParams{ll}.orig(1),yCoords' - fitParams{ll}.orig(2)];
        componentLines{ll} = componentLines{ll} * [cos(rotAngle(ll)),-sin(rotAngle(ll));sin(rotAngle(ll)),cos(rotAngle(ll))];
        componentLines{ll} = componentLines{ll} + repmat(mean(cc),[length(xCoords),1]);
        
        % 'componentLines{ll}' contains a string of coordinates, +/- 1 pixel apart,
        % identifying the location of the fish with subpixel accuracy.
        % The first and last points are not necessarily at the very
        % tips of the fish, so you can't use those to determine its
        % velocity. Just use these to measure shape.
        
        
        %             if invertedRotAngle==1
        %                 if (rotAngle(ll) + pi) > 2*pi
        %                     rotAngle(ll) = rotAngle(ll) - pi;
        %                 else
        %                     rotAngle(ll) = rotAngle(ll) + pi;
        %                 end
        %             end
        %
        
        
        
        
        
        
        
    end
    
    
    % match the new object to the old labels
    if (labelCounter == 0) % if this is the first frame, just give labels from one to N to all the N components
        objectLabels = 1:length(msr);
        previousObjectLabels = objectLabels;
        previousObjectCoordinates = msr.Center';
        labelCounter = labelCounter + length(msr);
    else
        
        %        previousObjectCoordinates
        %        previousObjectLabels
        currentObjectCoordinates = msr.Center';
        %        currentObjectCoordinates
        [objectLabels, nNewLabels] = match_coordinates(previousObjectCoordinates, currentObjectCoordinates, previousObjectLabels, thresholdDistance, labelCounter);
        labelCounter = labelCounter + nNewLabels;
        previousObjectLabels = objectLabels;
        previousObjectCoordinates = currentObjectCoordinates;
    end
    
    
    
end


% Display.
if (showFigures > 0)
    myFig = dipshow(imgDip);
    diptruesize(myFig, 50);
    hold on;
    myColours = 'rgbwyc';
    if ~isempty(msr)
        for ll = 1:length(msr)
            % figure(myFig);
            %            hold on;
            fitParams{ll}.validatedCurve = 0; % no curve is validated yet
            
            plot(componentLines{ll}(:,1),componentLines{ll}(:,2),myColours(mod(ll-1, 6)+1));
            plot(msr.Center(1,ll), msr.Center(2,ll), 'Marker', 'o', 'MarkerSize', 6, 'Color', 'k', 'MarkerFaceColor', [0.8 0.8 0.5]);
            text(componentLines{ll}(1,1) -20,componentLines{ll}(1,2)-20,num2str(objectLabels(ll)) ,'color','r');
            %            text(componentLines{ll}(end,1),componentLines{ll}(end,2),num2str(rotAngle(ll)*180/pi),'color','red')
            %            text(componentLines{ll}(end,1),componentLines{ll}(end,2)+20,num2str(fitParams{ll}.meanDeviation),'color','green')
            %            text(componentLines{ll}(end,1),componentLines{ll}(end,2)+40,num2str(fitParams{ll}.tailFinPosition),'color','green')
            %         text(componentLines{ll}(end,1),componentLines{ll}(end,2)+20,num2str(fitParams{ll}.phi*180/pi),'color','green')
            %         text(componentLines{ll}(end,1),componentLines{ll}(end,2)+40,num2str(fitParams{ll}.omega*180/pi),'color','blue')
        end
        
        %         objectContour = (bdilation(binaryImg) - binaryImg);
        %         overlay(myGrayImage, objectContour);
        %    selectedCoordinates = dipgetcoords(1);
        %    selectedObject = dip_array(labeledImage(selectedCoordinates(1), selectedCoordinates(2)))
        %
        %     figure,
        %     imgVisualCheck = dip_array(imgDip);
        %     componentsImage = dip_array(labeledImage);
        %     componentContours = bwmorph(componentsImage,'remove');
        %     colourMap = jet(length(msr))
        %     imgVisualCheck(find(componentContours ~= 0) = colourMap();
        %     hold on;
        %     image(dip_array(imgDip));
        
        
        validatedObjects = dipgetcoords(length(msr));
        validatedObjects(validatedObjects < 0) = inf;
        squaredEuclideanDistances = distm(validatedObjects(:,1:2), msr.Center');
        
        % Now find and validate the closest object to each left click
        [validatedDistance, validatedObject] = min(squaredEuclideanDistances, [], 2);
        for ll = 1:size(validatedObjects,1)
            if (validatedDistance(ll) < 10000)
                objectLabels(ll)
                fitParams{validatedObject(ll)}.validatedCurve = 1;
            end
        end
    end
    
    close(myFig);
    
end









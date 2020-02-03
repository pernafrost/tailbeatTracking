

clear all;

showFig = 0;


allData=[];
currElementFish1 = [];
currElementFish2 = [];
allAllFiles={};

allDirs = dir('input_images/');
for dd =  1:length(allDirs) % dd is the directory counter
    if allDirs(dd).name(1) ~= '.'
        allFiles = dir(fullfile('input_images', allDirs(dd).name, 'results*.mat'));
        allAllFiles{dd} = allFiles;
        for ff = 1:length(allFiles) % ff is the file counter
            inputFile = fullfile('input_images', allDirs(dd).name, allFiles(ff).name)
            clear results
            load(inputFile);
            
            
            Fs = 60; % hertz, sampling frequency
            
            
            maxNumObjects = 0;
            for jj =1:length(results) % jj is the frame counter
                labelValue = max(results{jj}.objectLabels);
                if labelValue > maxNumObjects
                    maxNumObjects = labelValue;
                end
            end
            
            %% select the good data and rearrange them "by fish"
            
            centreOfMass = cell(maxNumObjects, 1);
            fishMovement = cell(maxNumObjects, 1);
            tailFinPosition = cell(maxNumObjects, 1);
            fishFrames = cell(maxNumObjects, 1);
            for kk  = 1:maxNumObjects % kk is the fish counter
                for jj = 1:length(results) % jj is the frame counter
                    fishIndex = find (results{jj}.objectLabels == kk);
                    
                    if ~isempty(fishIndex)% results{frameCounter}.fitParams{1}.validatedCurve
                        
                        
                        % I want to exclude the fish that touch the edge of the image.
                        % This can be done either by looking at the minimum and maximum
                        % coordinates (take care that I do not measure the image size from the image
                        % but I assume that it is 1600x1200!) or focusing on the tip of
                        % the tail (even if the fish is not complete, it is fine as long
                        % as we can see the tip of the tail)
                        % the fish is valid if the tip of the tail is visible, as well
                        % as a proportion of its body
                        
                        tipOfTail = results{jj}.componentLines{fishIndex}(end,:);
                        isFishCutByImageFrame = (size(results{jj}.componentLines{fishIndex},1) < 140 || tipOfTail(1) <= 3 || tipOfTail(2) <= 3 || tipOfTail(1) >= 1600-4 || tipOfTail(2) >= 1200-4);
                        
                        if ~isFishCutByImageFrame
                            fishMovement{kk}(end+1, 1:7, 1) = results{jj}.fitParams{fishIndex}.lastSegment(1,:); % slope
                            fishMovement{kk}(end, 1:7, 2) = results{jj}.fitParams{fishIndex}.lastSegment(2,:); % intercept
                            tailFinPosition{kk}(end+1) = results{jj}.fitParams{fishIndex}.tailFinPosition; % position of the tip of the tail fin
                            centreOfMass{kk}(end+1,:) = results{jj}.msr.data{1}(:,fishIndex);
                            fishFrames{kk}(end+1) = jj;
                        end
                    end
                end
            end
            
            
            
            
            
            
            clear xLimits omega tailFinPositionNormalized speed xv yv movingDirection;
            for kk = 1:2
                
                
                medianAbsTailPosition = median(abs(tailFinPosition{kk}));
                % median(abs(sine wave)) = sqrt(2)/2 hence if I want my function to be
                % between -1 and 1 I have to multiply the median of absolute value by sqrt(2)
                tailFinPositionNormalized{kk} = tailFinPosition{kk}(:)/medianAbsTailPosition/sqrt(2);
                tailFinPositionDerivative = (tailFinPositionNormalized{kk}(3:end) - tailFinPositionNormalized{kk}(1:end-2))./(fishFrames{kk}(3:end)- fishFrames{kk}(1:end-2))';
                % for the first and last point I estimate the derivative only on one
                % side
                tailFinPositionDerivative = [diff(tailFinPositionNormalized{kk}(1:2))/diff(fishFrames{kk}(1:2)); ...
                    tailFinPositionDerivative; ...
                    diff(tailFinPositionNormalized{kk}(end-1:end))/diff(fishFrames{kk}(end-1:end));];
                
                xv = (centreOfMass{kk}(2:end,1) - centreOfMass{kk}(1:end-1,1)) ./(fishFrames{kk}(2:end)- fishFrames{kk}(1:end-1))';
                yv = (centreOfMass{kk}(2:end,2) - centreOfMass{kk}(1:end-1,2)) ./(fishFrames{kk}(2:end)- fishFrames{kk}(1:end-1))';
                [movingDirection{kk}, speed{kk}] = cart2pol(xv, yv);
                
                
                
                %     tailFinPositionDerivative = diff(tailFinPositionNormalized{kk})./diff(fishFrames{kk});
                omega{kk} = atan2(tailFinPositionNormalized{kk}, tailFinPositionDerivative); % the derivative of sine is the cosine
                
                % here test an alternative where omega is just the tailFinPositionNormalized
                % omega{kk} = tailFinPositionNormalized{kk};
                
            end
            
            
            markerColours = 'br';
            markerFaceColours = [0 1 1; 1 0.7 0.7];
            
            
            for jj =1:length(results)
                
                
                % display
                if showFig
                    % img = imread(results{jj}.inputImage);
                    
                    scrsz = get(0,'ScreenSize');
                    figure('Position', [1 1 scrsz(3)/2 scrsz(4)/2]);
                    hold on;
                    % figure('Position',[1 1 scrsz(3) scrsz(4)]);
                    % subplot(5,4,[1, 2, 5, 6]);
                    subplot(2,2,[1,3]);
                    % image(img);
                    axis image;
                    hold on;
                end
                
                if ~isempty(results{jj}.msr.data)
                    fish1Label = results{jj}.objectLabels(1); % label of first fish
                    diffCoordinates = diff(results{jj}.msr.data{1},[],2); % in the data each column is a fish
                    % the differnce in coordinates gives me the vector
                    % going from the first to the second fish; for the second fish I take
                    % opposite sign
                    
                    
                    
                    currElementFish1 = find(fishFrames{1} == jj, 1, 'first');
                    currElementFish2 = find(fishFrames{2} == jj, 1, 'first');
                    
                    
                   
                    
                    if showFig
                        plot(results{jj}.componentLines{1}(:,1), results{jj}.componentLines{1}(:,2));
                        plot(results{jj}.msr.data{1}(1,1), results{jj}.msr.data{1}(2,1), 'Marker', 'o' , 'Color', markerColours(mod(fish1Label-1,2) +1), 'MarkerSize', 7, 'MarkerFaceColor', markerFaceColours(mod(fish1Label-1,2) +1,:));
                        text(results{jj}.msr.data{1}(1,1)+20, results{jj}.msr.data{1}(2,1)+20, num2str(fish1Label));
                    end
                    
                    
                    
                    
                    if ~isempty(diffCoordinates)
                        fish2Label = results{jj}.objectLabels(2); % label of second fish
                        [angleFromHorizontal, r] = cart2pol([diffCoordinates(1); -diffCoordinates(1)], [diffCoordinates(2); - diffCoordinates(2)]);
                        % the direction of movement of a fish is
                        fishOwnDirection = pi - results{jj}.rotAngle;
                        
                        % the distance between two neighbours is
                        theta = mod(pi + angleFromHorizontal - fishOwnDirection, 2*pi) - pi;
                        phi = [mod(pi + diff(fishOwnDirection), 2*pi) - pi; mod(pi - diff(fishOwnDirection), 2*pi) - pi];
                        
                        oneLineOfData(1) = dd;
                        oneLineOfData(2) = ff;
                        oneLineOfData(3) = jj;
                        oneLineOfData(4) = fish1Label;
                        oneLineOfData(5) = fishOwnDirection(1);
                        oneLineOfData(6) = theta(1);
                        oneLineOfData(7) = phi(1);
                        oneLineOfData(8) = r(1);
                        if isempty(currElementFish1) || (currElementFish1 == 1)
                            oneLineOfData(9) = NaN;
                            oneLineOfData(10) = NaN;
                        else
                            oneLineOfData(9) = movingDirection{1}(currElementFish1-1);
                            oneLineOfData(10) = speed{1}(currElementFish1-1);
                        end
                        if isempty (currElementFish1)
                            oneLineOfData(11) = NaN;
                        else
                            oneLineOfData(11) = omega{1}(currElementFish1);
                        end
                        if isempty (currElementFish2)
                            oneLineOfData(12) = NaN;
                        else
                            oneLineOfData(12) = omega{2}(currElementFish2);
                        end
                        
                        
                        allData = [allData; oneLineOfData];
                        
                        
                        oneLineOfDataOld = oneLineOfData;
                        
                        oneLineOfData(4) = fish2Label;
                        oneLineOfData(5) = fishOwnDirection(2);
                        oneLineOfData(6) = theta(2);
                        oneLineOfData(7) = phi(2);
                        oneLineOfData(8) = r(2);
                        if isempty(currElementFish2) || (currElementFish2 == 1)
                            oneLineOfData(9) = NaN;
                            oneLineOfData(10) = NaN;
                        else
                            oneLineOfData(9) = movingDirection{2}(currElementFish2-1);
                            oneLineOfData(10) = speed{2}(currElementFish2-1);
                        end
                        if isempty (currElementFish2)
                            oneLineOfData(11) = NaN;
                        else
                            oneLineOfData(11) = omega{2}(currElementFish2);
                        end
                        if isempty (currElementFish1)
                            oneLineOfData(12) = NaN;
                        else
                            oneLineOfData(12) = omega{1}(currElementFish1);
                        end
                        
                        
                        
                        allData = [allData; oneLineOfData];
                        
                        
                        
                        if showFig
                            % plot(results{jj}.msr.data{1}(1,:), results{jj}.msr.data{1}(2,:), 'Marker', 'o');
                            % text(results{jj}.msr.data{1}(1,:)+20, results{jj}.msr.data{1}(2,:)+20, {'1', '2'});
                            plot(results{jj}.msr.data{1}(1,2), results{jj}.msr.data{1}(2,2), 'Color', markerColours(mod(fish2Label-1,2) +1), 'MarkerSize', 7, 'MarkerFaceColor', markerFaceColours(mod(fish2Label-1,2) +1,:));
                            text(results{jj}.msr.data{1}(1,2)+20, results{jj}.msr.data{1}(2,2)+20, num2str(fish2Label));
                            arrowToNeighbour = quiver(results{jj}.msr.data{1}(1,:)', results{jj}.msr.data{1}(2,:)', 100*cos(angleFromHorizontal), 100*sin(angleFromHorizontal), 0);
                            set(arrowToNeighbour, 'LineWidth',2);
                            arrowToOwnDirection = quiver(results{jj}.msr.data{1}(1,:)', results{jj}.msr.data{1}(2,:)', 100*cos(fishOwnDirection), 100*sin(fishOwnDirection), 0);
                            plot(results{jj}.componentLines{2}(:,1), results{jj}.componentLines{2}(:,2));
%                             oneLineOfDataOld
%                             oneLineOfData
%                             pause (0.4);
                        end
                        
                        
                        
                    end
                    %
                    % subplot(5,4,[9, 10, 13, 14]);
                    % polar(3,2);
                    % % text('rotated');
                    
                    
                end
                
                
                
                

                
                
                if showFig
                    
                    subplot(2,2,2);
                    plot(fishFrames{1}/Fs, tailFinPositionNormalized{1}, 'c', 'LineWidth', 1); % , 'LineStyle', 'none', 'Marker', '.');
                    
                    if ~isempty(currElementFish1),
                        hold on;
                        plot(fishFrames{1}(currElementFish1)/Fs, tailFinPositionNormalized{1}(currElementFish1), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'b');
                    end
                    set(gca, 'XLim', [0, max(fishFrames{1}(end), fishFrames{2}(end))/Fs]);
                    
                    subplot(2,2,4);
                    plot(fishFrames{2}/Fs, tailFinPositionNormalized{2}, 'Color', [1 0.7 0.7], 'LineWidth', 1); %, 'LineStyle', 'none', 'Marker', '.');
                    
                    if ~isempty(currElementFish2),
                        hold on;
                        plot(fishFrames{2}(currElementFish2)/Fs, tailFinPositionNormalized{2}(currElementFish2), 'LineStyle', 'none', 'Marker', 'o', 'Color', 'r');
                    end
                    set(gca, 'XLim', [0, max(fishFrames{1}(end), fishFrames{2}(end))/Fs]);
                    
                    figure('Position', [scrsz(3)*5/9 1 scrsz(3)/3 scrsz(4)/2]);
                    
                    
                    
                    if ~isempty(currElementFish1)
                        polarPlot1 = polar(omega{1}(currElementFish1), 1);
                        set(polarPlot1, 'Marker', 'o', 'Color', 'b', 'MarkerSize', 12, 'MarkerFaceColor', 'c');
                        hold on;
                    end
                    
                    if ~isempty(currElementFish2)
                        polarPlot2 = polar(omega{2}(currElementFish2), 1);
                        set(polarPlot2, 'Marker', 'o', 'Color', 'r', 'MarkerSize', 12, 'MarkerFaceColor', [1 0.7 0.7]);
                        hold on;
                    end
                    
                    pause(0.4);
                    close all;
                    
                end
            end
            
            
            % dir number, video number, frame number, fish number,
            % speed, acceleration, turn, r, theta, phi, speed
            % neighbour, omega focal, omega neighbour
            
            % also save dir number, video number, fish number, size, peak
            % frequency, second peak frequency, amplitude, second peak
            % amplitude in a different array
            
            
        end
    end
end



save body_movement_analysed_data.mat allData allDirs allAllFiles;


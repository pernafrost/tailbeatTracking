

clear all;

allData = [];

allDirs = dir('input_images');
for dd = 1:length(allDirs) % dd is the directory counter
    if allDirs(dd).name(1) ~= '.'
        allFiles = dir(fullfile('input_images', allDirs(dd).name, 'results*.mat'));
        for ff = 1:length(allFiles) % ff is the file counter
            clear results
            load(fullfile('input_images', allDirs(dd).name, allFiles(ff).name))
            fullfile(allDirs(dd).name, allFiles(ff).name)
            % load images/T1/results_11-32-17.726-frames_133453-133522.mat
            
            Fs = 60; % hertz, sampling frequency
            
            
            maxNumObjects = 0;
            for jj =1:length(results) % jj is the frame counter
                labelValue = max(results{jj}.objectLabels);
                if labelValue > maxNumObjects
                    maxNumObjects = labelValue;
                end
            end
            
            %% select the good data and rearrange them "by fish"
            
            fishMovement = cell(maxNumObjects, 1);
            tailFinPosition = cell(maxNumObjects, 1);
            fishFrames = cell(maxNumObjects, 1);
            fishSize = cell(maxNumObjects, 1);
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
                        
                        
                        %             fishMinima = strmatch('Minimum', results{jj}.msr.names, 'exact');
                        %             fishMaxima = strmatch('Maximum', results{jj}.msr.names, 'exact');
                        %             fishXMin = results{jj}.msr.data{fishMinima}(1, fishIndex);
                        %             fishYMin = results{jj}.msr.data{fishMinima}(2, fishIndex);
                        %             fishXMax = results{jj}.msr.data{fishMaxima}(1, fishIndex);
                        %             fishYMax = results{jj}.msr.data{fishMaxima}(2, fishIndex);
                        %
                        %             % check if the fish is touching the edge of the image (in which case I do not consider it for some of the analyses)
                        %            isFishCutByImageFrame = (fishXMin == 0 || fishYMin == 0 || fishXMax == 1600 || fishYMax == 1200);
                        %
                        
                        
                        
                        
                        % the fish is valid if the tip of the tail is visible, as well
                        % as a proportion of its body
                        %            jj
                        %            fishIndex
                        tipOfTail = results{jj}.componentLines{fishIndex}(end,:);
                        isFishCutByImageFrame = (size(results{jj}.componentLines{fishIndex},1) < 140 || tipOfTail(1) <= 3 || tipOfTail(2) <= 3 || tipOfTail(1) >= 1600-4 || tipOfTail(2) >= 1200-4);
                        
                        % check that I am doing things properly
                        
                        %            img = imread(results{jj}.inputImage);
                        %            image(img); hold on; plot(results{jj}.componentLines{fishIndex}(:,1), results{jj}.componentLines{fishIndex}(:,2));
                        %            pause; close all;
                        
                        
                        
                        if ~isFishCutByImageFrame
                            
                            fishMovement{kk}(end+1, 1:7, 1) = results{jj}.fitParams{fishIndex}.lastSegment(1,:); % slope
                            fishMovement{kk}(end, 1:7, 2) = results{jj}.fitParams{fishIndex}.lastSegment(2,:); % intercept
                            tailFinPosition{kk}(end+1) = results{jj}.fitParams{fishIndex}.tailFinPosition; % position of the tip of the tail fin
                            % fishSize{kk}(end+1) = results{jj}.msr.data{strmatch('Size', results{jj}.msr.names)}(fishIndex);
                            fishSize{kk}(end+1) = results{jj}.msr(fishIndex).Size;

                            fishFrames{kk}(end+1) = jj;
                        end
                    end
                end
            end
            
            
            
            
            
            
            
            
            
            %% Interpolate data if they are missing (there shouldn't be missing data).
            
            for kk = 1:2
                if ~isempty(fishFrames{kk})
                tailFinPositionInterp{kk} = interp1(fishFrames{kk}, tailFinPosition{kk}, min(fishFrames{kk}):max(fishFrames{kk}), 'spline');
            
                end
            end
            
            
            %% try fourier transform
            
            
            
            colours = 'rgbk';
            
            clear peakFrequency;
            scrsz = get(0,'ScreenSize');
            figure('Position',[1 scrsz(4)/10 scrsz(3)*3/4 scrsz(4)*3/4]);
            for kk = 1:2
                
                %    for ii = 1:4
                ii = 1;
                % myData = fishMovement{jj}(isfinite(fishMovement{jj}(:,ii,1)),ii,1); % put to zero the values that are not finite
                myData = tailFinPositionInterp{kk};
                NFFT = 2^nextpow2(length(myData)); % Next power of 2 from length of y
                myDataFreq = fft(myData,NFFT)/length(myData);
                f = Fs/2*linspace(0,1,NFFT/2+1);
                
                
                % Plot single-sided amplitude spectrum.
                subplot(2,3, 1 + (kk-1)*3);
                plot(f,2*abs(myDataFreq(1:floor(NFFT/2+1))), colours(mod(ii-1,4)+1));
                title('Single-Sided Amplitude Spectrum of fish tail movement')
                xlabel('Frequency (Hz)')
                ylabel('|myDataFreq(f)|');
                hold on;
                
                [peakAmplitude(kk),I] = max(2*abs(myDataFreq(1:floor(NFFT/2+1))));
                peakFrequency(kk) = f(I);
                
                % second peak
                myDataFreq(I) = 0;
                [peakAmplitude(kk+2),I] = max(2*abs(myDataFreq(1:floor(NFFT/2+1))));
                peakFrequency(kk+2) = f(I);
                
                plot(peakFrequency(kk), peakAmplitude(kk), 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 12 , 'Color', 'k');
                
                plot(peakFrequency(kk+2), peakAmplitude(kk+2), 'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 12, 'Color', 'b');
                %    end
            end
            
            
            
            % prctile(abs(tailFinPosition{kk}),90)
            % median(abs(tailFinPosition{kk}))
            
            
            
            %% plot the data
            
            
            clear xLimits omega tailFinPositionNormalized;
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
                
                
                %     tailFinPositionDerivative = diff(tailFinPositionNormalized{kk})./diff(fishFrames{kk});
                omega{kk} = atan2(tailFinPositionNormalized{kk}, tailFinPositionDerivative);
                
                subplot(2,3, [2 + (kk-1)*3, 3 + (kk-1)*3]);
                
                hold on;
                % plot(fishFrames{jj}(:), fishMovement{jj}(:,1,1), 'r');
                % plot(fishFrames{jj}(:), fishMovement{jj}(:,2,1), 'g');
                % plot(fishFrames{jj}(:), fishMovement{jj}(:,3,1), 'b');
                % plot(fishFrames{jj}(:), fishMovement{jj}(:,4,1), 'c');
                % plot(fishFrames{jj}(:), fishMovement{jj}(:,5,1), 'k');
                % plot(fishFrames{jj}(:), fishMovement{jj}(:,6,1), 'y');
                plot(fishFrames{kk}(:)/Fs, tailFinPositionNormalized{kk}, 'c', 'LineWidth', 2);
                %     plot((0:0.1:60)/Fs, peakAmplitude(kk)*sin((0:0.1:60)/Fs*peakFrequency(kk)*2*pi), 'r', 'LineWidth', 1.4);
                %     plot(fishFrames{kk}(:)/Fs, peakAmplitude(kk)*sin(fishFrames{kk}(:)/Fs*peakFrequency(kk)*2*pi), 'r', 'LineWidth', 1.4);
                xLimits(kk,:) = get(gca, 'XLim');
                %    boxplot(tailFinPosition{kk}, 'positions', 0.2, 'widths', 0.1);
            end
            
            
            for kk = 1:2
                subplot(2,3,[2 + (kk-1)*3, 3 + (kk-1)*3]);
                set(gca, 'XLim', [min(xLimits(:,1),[],1), max(xLimits(:,2),[],1)]);
                set(gca, 'XTick', 0:0.2:1);
                plot([min(xLimits(:,1),[],1), max(xLimits(:,2),[],1)], [0,0], 'k');
                xlabel('Time (s)');
                ylabel('tail displacement (pixels)');
            end
            
            
            for kk = 1:2
                meanSize = mean(fishSize{kk});
                stdSize = std(fishSize{kk});
                nFrames = length(fishSize{kk});
                
                oneLineOfData = [dd, kk, meanSize, stdSize, nFrames, peakFrequency(kk), peakFrequency(kk+2), peakAmplitude(kk), peakAmplitude(kk+2)];
                allData = [allData; oneLineOfData];
            end
                
            pause(2);
            close all;
            
        end
    end
end


save data_on_frequency.mat allData;


load data_on_frequency.mat
dd = 7;
aaa = find(allData(:,1) == dd); %  find data in the same folder (same pair of fish)
figure, plot(allData(aaa,3), allData(aaa,6), '.')

figure, plot(allData(:,3), allData(:,6), '.')
xlabel('blob size'); ylabel('peak frequency');


figure, plot(allData(1:2:end,6), allData(2:2:end,6), '.');
xlabel('freq. fish 1'); ylabel('freq. fish 2');
axis square; % axis equal;
set(gca,'XLim', [0, 15], 'YLim', [0, 15]);

% allData = allData(allData(:,6) > 5, :);





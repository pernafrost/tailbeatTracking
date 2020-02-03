

clear all;


load body_movement_analysed_data.mat

                   % allData(1) = dd;
                   % allData(2) = ff;
                   % allData(3) = jj;
                   % allData(4) = fish1Label;
                   % allData(5) = fishOwnDirection(1);
                   % allData(6) = theta(1);
                   % allData(7) = phi(1);
                   % allData(8) = r(1);
                   % allData(9) = movingDirection{1}(currElementFish1-1);
                   % allData(10) = speed{1}(currElementFish1-1);
                   % allData(11) = omega{1}(currElementFish1);
                   % allData(12) = omega{2}(currElementFish2);


deltaOmega = mod(allData(:,12)- allData(:,11) + pi, 2*pi)-pi;

% make a colourful plot
distanceBins = 0:50:1000;
thetaBins = linspace(-pi, pi, 13);

[nInD, valInD] = histc(allData(:,8), distanceBins);
[nInTheta, valInTheta] = histc(allData(:,6), thetaBins);



counts = zeros(length(distanceBins), length(thetaBins));
density = zeros(length(distanceBins), length(thetaBins));
meanOmega = zeros(length(distanceBins), length(thetaBins));
circularStandardDeviationOmega = zeros(length(distanceBins), length(thetaBins));

for countD=1:length(distanceBins)
    for countTheta=1:length(thetaBins)
        valDTheta = intersect(find(valInD==countD), find(valInTheta==countTheta));
        if countD < length(distanceBins)
            area = pi/6 * (distanceBins(countD+1)^2 - distanceBins(countD)^2) / 2;
        else % for the last value area is not important because the value is nan (?)
            area = 1;
        end
        counts(countD, countTheta) = length(valDTheta);
        density(countD, countTheta) = counts(countD, countTheta)./area;
        
        xProjOmega = nanmean(cos(deltaOmega(valDTheta)));
        % xProjOmega = nanmean(cos(deltaOmega(valDTheta)).*sign(allData(valDTheta, 6)));
        yProjOmega = nanmean(sin(deltaOmega(valDTheta)));
        r = sqrt(xProjOmega^2 + yProjOmega^2);
        meanOmega(countD,countTheta) = atan2(yProjOmega, xProjOmega);
        circularStandardDeviationOmega(countD,countTheta) = sqrt(-2*log(r));
    end
end



X = distanceBins'*cos(pi/2 +thetaBins);
Y = distanceBins'*sin(pi/2 +thetaBins);

useFontSize = 18;

valueForColourScale = max(max(density));
figure,
set(gcf, 'Position', [1 1 1000 640]);
pcolor(X,Y,density);
caxis([0, valueForColourScale]);
axis equal xy tight;
box off;
colormap(jet);
set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'TickDir', 'out', 'XTick', -1000:100:1000, 'YTick', -1000:100:1000, ...
    'XLim', [-1200 1200], 'YLim', [-1200 1200], 'XTickLabel', -1000:100:1000, 'YTickLabel', -1000:100:1000, 'FontName', 'Arial');
text(0,max(distanceBins)*11/10,'\vartheta = 0', 'FontName', 'Arial', 'FontSize', useFontSize, 'VerticalAlignment','middle', 'HorizontalAlignment','center');
text(max(distanceBins)*11/10,0,'\pi/2', 'FontName', 'Arial', 'FontSize', useFontSize, 'VerticalAlignment','middle', 'HorizontalAlignment','center');
text(-max(distanceBins)*11/10,0,'-\pi/2', 'FontName', 'Arial', 'FontSize', useFontSize, 'VerticalAlignment','middle', 'HorizontalAlignment','center');
text(0,-max(distanceBins)*11/10,'\pi', 'FontName', 'Arial', 'FontSize', useFontSize, 'VerticalAlignment','middle', 'HorizontalAlignment','center');
xlabel('distance from neighbour r (pixels)'); ylabel('distance from neighbour r (pixels)');

text(-1000, 1000, 'A', 'FontName', 'Arial', 'FontSize', 32);

% title('mean acceleration as a function of distance and angle of all neighbours');
colorbarHandle = colorbar;
set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial');%  , 'YTick', 0:0.1:1.2, 'YTickLabel', 0:0.1:1.2);
ylabel(colorbarHandle, 'density (counts/pixel^2)', 'FontSize', useFontSize);







figure,
set(gcf, 'Position', [1 1 1000 640]);
pcolor(X,Y,meanOmega);
caxis([-pi, pi]);
axis equal xy tight;
box off;
colormap(jet);
set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'TickDir', 'out', 'XTick', -1000:100:1000, 'YTick', -1000:100:1000, ...
    'XLim', [-1200 1200], 'YLim', [-1200 1200], 'XTickLabel', -1000:100:1000, 'YTickLabel', -1000:100:1000, 'FontName', 'Arial');
text(0,max(distanceBins)*11/10,'\vartheta = 0', 'FontName', 'Arial', 'FontSize', useFontSize, 'VerticalAlignment','middle', 'HorizontalAlignment','center');
text(max(distanceBins)*11/10,0,'\pi/2', 'FontName', 'Arial', 'FontSize', useFontSize, 'VerticalAlignment','middle', 'HorizontalAlignment','center');
text(-max(distanceBins)*11/10,0,'-\pi/2', 'FontName', 'Arial', 'FontSize', useFontSize, 'VerticalAlignment','middle', 'HorizontalAlignment','center');
text(0,-max(distanceBins)*11/10,'\pi', 'FontName', 'Arial', 'FontSize', useFontSize, 'VerticalAlignment','middle', 'HorizontalAlignment','center');
xlabel('distance from neighbour r (pixels)'); ylabel('distance from neighbour r (pixels)');

text(-1000, 1000, 'A', 'FontName', 'Arial', 'FontSize', 32);

% title('mean acceleration as a function of distance and angle of all neighbours');
colorbarHandle = colorbar;
set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial'); % , 'YTick', 0:0.1:1.2, 'YTickLabel', 0:0.1:1.2);
ylabel(colorbarHandle, 'meanOmega (counts/pixel^2)', 'FontSize', useFontSize);












%%% IF omega is an angle



% make a colourful plot
xDistanceBins = -500:50:500;
yDistanceBins = -500:50:500;

% [nInXD, valInXD] = histc(allData(:,8).*abs(sin(allData(:,6))), xDistanceBins);
[nInXD, valInXD] = histc(allData(:,8).*sin(allData(:,6)), xDistanceBins);
[nInYD, valInYD] = histc(allData(:,8).*cos(allData(:,6)), yDistanceBins);

% only take into account fish that are more or less aligned to each other
alignedFish = find(abs(allData(:,7)) < pi/8);


counts = zeros(length(yDistanceBins), length(xDistanceBins));
meanOmega = zeros(length(yDistanceBins), length(xDistanceBins));
circularStandardDeviationOmega = zeros(length(yDistanceBins), length(xDistanceBins));

for countXD=1:length(xDistanceBins)
    for countYD=1:length(yDistanceBins)
        % valXY = intersect(find(valInXD==countXD), find(valInYD==countYD));
        valXY = intersect(alignedFish, intersect(find(valInXD==countXD), find(valInYD==countYD)));
        
        counts(countYD, countXD) = length(valXY);
        
        xProjOmega = nanmean(cos(deltaOmega(valXY)));
        yProjOmega = nanmean(sin(deltaOmega(valXY)));
        r = sqrt(xProjOmega^2 + yProjOmega^2);
        meanOmega(countYD,countXD) = atan2(yProjOmega, xProjOmega);
        circularStandardDeviationOmega(countYD,countXD) = sqrt(-2*log(r));
    end
end


meanOmega(counts<20) = 0;

figure,
set(gcf, 'Position', [1 1 1000 640]);
imagesc(meanOmega(1:end-1, 1:end-1));
caxis([-pi, pi]);
% caxis([-1, 1]);
axis equal xy tight;
box off;
colormap(hsv);
 set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'XTick', 0.5:2:length(xDistanceBins), 'YTick', 0.5:2:length(yDistanceBins),...
      'XTickLabel', xDistanceBins(1:2:end), 'YTickLabel', yDistanceBins(1:2:end), 'TickDir', 'out', ...
        'FontName', 'symbol');
    axis xy equal tight;
    xlabel('X distance (pixels)', 'FontName', 'arial'); ylabel('Y distance (pixels)', 'FontName', 'arial');

    % title('mean acceleration as a function of distance and angle of all neighbours');
colorbarHandle = colorbar;
set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial'); % , 'YTick', 0:0.1:1.2, 'YTickLabel', 0:0.1:1.2);
ylabel(colorbarHandle, 'mean delta omega', 'FontSize', useFontSize);



    
    
figure,
set(gcf, 'Position', [1 1 1000 640]);
imagesc(counts(1:end-1, 1:end-1));
caxis([0, max(max(counts))]);
axis equal xy tight;
box off;
colormap(jet);
 set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'XTick', 0.5:2:length(xDistanceBins), 'YTick', 0.5:2:length(yDistanceBins),...
      'XTickLabel', xDistanceBins(1:2:end), 'YTickLabel', yDistanceBins(1:2:end), 'TickDir', 'out', ...
        'FontName', 'symbol');
    axis xy equal tight;
    xlabel('X distance (pixels)', 'FontName', 'arial'); ylabel('Y distance (pixels)', 'FontName', 'arial');

    % title('mean acceleration as a function of distance and angle of all neighbours');
colorbarHandle = colorbar;
set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial'); % , 'YTick', 0:0.1:1.2, 'YTickLabel', 0:0.1:1.2);
ylabel(colorbarHandle, 'number of counts', 'FontSize', useFontSize);



    
    
    

% 
% 
% 
% %%% IF omega is a position
% 
% 
% 
% % make a colourful plot
% xDistanceBins = -500:50:500;
% yDistanceBins = -500:50:500;
% 
% [nInXD, valInXD] = histc(allData(:,8).*abs(cos(allData(:,6))), xDistanceBins);
% % [nInXD, valInXD] = histc(allData(:,8).*cos(allData(:,6)), xDistanceBins);
% [nInYD, valInYD] = histc(allData(:,8).*sin(allData(:,6)), yDistanceBins);
% 
% 
% 
% counts = zeros(length(yDistanceBins), length(xDistanceBins));
% meanOmega = zeros(length(yDistanceBins), length(xDistanceBins));
% standardDeviationOmega = zeros(length(yDistanceBins), length(xDistanceBins));
% 
% for countXD=1:length(xDistanceBins)
%     for countYD=1:length(yDistanceBins)
%         valXY = intersect(find(valInXD==countXD), find(valInYD==countYD));
%         
%         counts(countYD, countXD) = length(valXY);
%         
%         % xProjOmega = nanmean(cos(deltaOmega(valXY)));
%         % yProjOmega = nanmean(sin(deltaOmega(valXY)));
%         % r = sqrt(xProjOmega^2 + yProjOmega^2);
%         meanOmega(countYD,countXD) = nanmean(deltaOmega(valXY).*sign(cos(allData(valXY,6))));% atan2(yProjOmega, xProjOmega);
%         standardDeviationOmega(countYD,countXD) = nanstd(deltaOmega(valXY));
%     end
% end
% 
% 
% meanOmega(counts<10) = 0;
% 
% figure,
% set(gcf, 'Position', [1 1 1000 640]);
% imagesc(meanOmega(1:end-1, 1:end-1));
% % caxis([-pi, pi]);
% caxis([-1, 1]);
% axis equal xy tight;
% box off;
% colormap(jet);
%  set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'XTick', 0.5:2:length(xDistanceBins), 'YTick', 0.5:2:length(yDistanceBins),...
%       'XTickLabel', xDistanceBins(1:2:end), 'YTickLabel', yDistanceBins(1:2:end), 'TickDir', 'out', ...
%         'FontName', 'symbol');
%     axis xy equal tight;
%     xlabel('X distance (pixels)', 'FontName', 'arial'); ylabel('Y distance (pixels)', 'FontName', 'arial');
% 
%     % title('mean acceleration as a function of distance and angle of all neighbours');
% colorbarHandle = colorbar;
% set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial'); % , 'YTick', 0:0.1:1.2, 'YTickLabel', 0:0.1:1.2);
% ylabel(colorbarHandle, 'mean delta omega', 'FontSize', useFontSize);
% 
% 
% 
%     
%     
% figure,
% set(gcf, 'Position', [1 1 1000 640]);
% imagesc(counts(1:end-1, 1:end-1));
% caxis([0, max(max(counts))]);
% axis equal xy tight;
% box off;
% colormap(jet);
%  set(gca, 'LineWidth', 2, 'FontSize', useFontSize, 'XTick', 0.5:2:length(xDistanceBins), 'YTick', 0.5:2:length(yDistanceBins),...
%       'XTickLabel', xDistanceBins(1:2:end), 'YTickLabel', yDistanceBins(1:2:end), 'TickDir', 'out', ...
%         'FontName', 'symbol');
%     axis xy equal tight;
%     xlabel('X distance (pixels)', 'FontName', 'arial'); ylabel('Y distance (pixels)', 'FontName', 'arial');
% 
%     % title('mean acceleration as a function of distance and angle of all neighbours');
% colorbarHandle = colorbar;
% set(colorbarHandle, 'FontSize', useFontSize, 'FontName', 'Arial'); % , 'YTick', 0:0.1:1.2, 'YTickLabel', 0:0.1:1.2);
% ylabel(colorbarHandle, 'number of counts', 'FontSize', useFontSize);
% 
% 
% 
%     
%     
%     
% 

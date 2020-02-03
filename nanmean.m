function [myMeanOmittingNaN] = nanmean(myInputVector)
%The function nanmean comes with the Statistics and Machine Learning
%Toolbox, but as I do not have the toolbox I implement a version here.

myMeanOmittingNaN = mean(myInputVector, 'omitnan');

end


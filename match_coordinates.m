function [proposedLabels, nNewLabels] = match_coordinates(previousObjectCoordinates, currentObjectCoordinates, previousObjectLabels, thresholdDistance, labelCounter)



D = distm(previousObjectCoordinates,currentObjectCoordinates);
assignedNewCoordinates = zeros(size(D, 2), 1);
assignedPreviousCoordinates = zeros(size(D, 1), 1);
proposedLabels = zeros(size(D, 2), 1);
Dpositions = 1:numel(D);
[rowIndices, columnIndices] = ind2sub(size(D), Dpositions);
D1(:,1) = D(Dpositions);
D1(:,2) = rowIndices; % each row is one of the previous coordinates
D1(:,3) = columnIndices; % each column is one of the current coordinates
D1 = sortrows(D1,1);

countMatch=0;
while nnz(assignedNewCoordinates) <= size(D,1) && (countMatch < size(D1,1))
    countMatch=countMatch+1;
    if D1(countMatch,1) > thresholdDistance^2
        break;
    end
    if ((assignedPreviousCoordinates(D1(countMatch,2))==0) && (assignedNewCoordinates(D1(countMatch,3))==0))
        assignedPreviousCoordinates(D1(countMatch,2)) = D1(countMatch,3);
        assignedNewCoordinates(D1(countMatch,3)) = D1(countMatch,2);
    end
end

nNewLabels = length(find(assignedNewCoordinates == 0));
proposedLabels(assignedNewCoordinates ~= 0) = previousObjectLabels(assignedNewCoordinates(assignedNewCoordinates ~= 0));
proposedLabels(assignedNewCoordinates == 0) = (1:nNewLabels) + labelCounter;



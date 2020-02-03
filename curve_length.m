function length = curve_length(xCoords,yCoords)

% syms x
% % Using the calculus formula
% % S=int(sqrt(1+dy/dx^2),a,b)
% % finding dy/dx
% poly_dif=diff(yCoords,x,1);
% % applying the formula
% integrand=sqrt(1+poly_dif^2);
% length=double(int(integrand,x,min(xCoords),max(xCoords)));
% 

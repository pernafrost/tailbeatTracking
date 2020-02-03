function out = fitting_my_function(params,data)
% scale x so that it goes from zero to one and from the tip of the tail
% from the center of the fish
% t = 2*pi*(x - xmin) / (xmax - xmin)

% scale y so that it is centered around zero and normalized in the same way
% as x with the (half) length of the fish
% q =  2*(y - ymean) / (ymax - ymin)

% Then fit the following function to the fish line
% q = cos({({omega t^{gamma}+phi})})*exp({-beta t})

% y=cos({({5x^{0.9}+1})})|_cdot_exp({-1.5x})

% x and y (t and q) are known
% find omega, gamma, beta and phi. phi is the only parameter which
% eventually changes once the others are known.


t = pi*(data(:,1) - min(data(:,1))) / (max(data(:,1)) - min(data(:,1))); % I multiply by pi and not 2*pi because x goes from the tip of the tail to the centre of the fish
q = (data(:,2) - mean(data(:,2))) / (max(data(:,1)) - min(data(:,1))); % amplitude normalized by fish length
phi = params(1);
omega = params(2);
A = params(3);
gamma = params(4);
beta = params(5);

% functionOft = A*cos(omega*t.^gamma + phi).*exp(-beta*t);
functionOft = A*cos(omega*t.^gamma + phi).*exp(-beta*t);


% minimum squared distance between observed and predicted
out = sum((q - functionOft).^2);


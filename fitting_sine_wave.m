function out = fitting_sine_wave(params,data)

% fit something like:
% y = A*sin(omega*t + phi);

t = data(1,:);
y = data(2,:);

A = params(1);
omega = params(2);
phi = params(3);
% DC = params(4);
% k = params(5); % linear trend



functionOft = A*sin(omega*t + phi); % + DC + k*t;


% minimum squared distance between observed and predicted
out = sum((y - functionOft).^2);


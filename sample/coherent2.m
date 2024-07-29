clear all
private_addpath('Advanpix/');

dim = 100;
domain = mp('[-pi pi;-pi pi]');
%domain = [0 pi;0 pi];

k = 1;
tau = 1; %mp('0.1');
basis = 'p';

T = @(x) x.^2/2;
V = @(x) cos(x);
dT = @(x) x;
dV = @(x) -sin(x);

[sU, state] = SplitUnitary(dim, domain, basis);

s.hsmrep()

[x,y,z] = s.hsmrep('ismp', true);
contour(ax, x, y, log10(z), 10 , 'LineColor', 'none', 'Fill','on');
clear all
e
%private_addpath('Advanpix');

%mp.Digits(200);

%sprintf("splitting.dat")
k = 1/2; %mp('1/2');
basis = 'p';
T = @(x) x.^2/2;
%v = @(x) atan( k * sin(2*pi*x) ./ (1 + k*cos(2*pi*x)))/pi;
%V = @(x) arrayfun(@(x) integral(v, 0, x), x);
V = @(x) cos(2*pi*x);
%domain = %mp('[-1 1;-1/2 1/2]');
domain = [-1/2 3/2; -1/2 1/2];

siorder = 1;
tau = 1;
dim = 100;
dirname = sprintf("test_N%d",dim);

sH = SplitHamiltonian(dim, domain, basis);
sU = SplitUnitary(dim, domain, basis);
hbar = sU.hbar;

mT = sH.matT(T);
mV = sH.matV(V);
tic 
disp("eig Hamiltonian")
Phi = mT + mV;

[pevecs, pevalsmat] = eig(Phi);
[pevals, sindex] = sort(real(diag(pevalsmat)), 'descend');
pevecs = pevecs(:, sindex);
pstates= eigs2states(sH, pevecs, pevals);
%utils.saveeigs(sH, pevals, pstates, 'basis', 'q', 'header', 'ham', 'savedir', dirname);
pslp = abs(pevals(2) - pevals(1));

toc 
disp("const matU")
tic
matU = sU.SImatrix(T, V, 'tau', tau, 'order', siorder);
toc

tic
disp("eig. matU")
[uevecs, uevalsmat] = eig(matU);
toc

tic
disp("sort matU")
%sindex = sU.sortbynorm(uevecs, pevecs);
%uevecs = uevecs(:, sindex);
%uevals = diag( uevals(sindex, sindex) );
uevals = diag( uevalsmat );
toc

ustates = eigs2states(sU, uevecs, uevals);
utils.saveeigs(sU, uevals, ustates, 'basis', 'q', 'header', 'uni', 'savedir', dirname);

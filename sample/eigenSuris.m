clear all
addpath('/Users/hanada/Dropbox/Packages/qtmaplab/');
private_addpath('Advanpix');

dim = 50; %mp.Digits(100);

domain = [0 1;-1/2 1/2]
domain = double(domain);

basis = 'p';
sH = SplitHamiltonian(dim, domain, basis);
sU = SplitUnitary(dim, domain, basis);

k = 1/2;
T = @(x) x.^2/2;
v = @(x) atan( k * sin(2*pi*x) ./ (1 + k*cos(2*pi*x)));
V = @(x) arrayfun(@(x) integral(v, 0, x, 'RelTol',0,'AbsTol',1e-16), x);

mCosP = sH.matT(@(x) cos(2*pi*x));
mSinP = sH.matT(@(x) sin(2*pi*x));
mCosQ = sH.matV(@(x) cos(2*pi*x));
mSinQ = sH.matV(@(x) sin(2*pi*x));
WeylQntz = @(x,y) (x * y + y * x)/2

siorder = 1;
tau = 1;
tau = mp('1');

tic 
disp("eig Hamiltonian")
Phi = mCosP + k * ( mCosQ + WeylQntz(mCosQ, mCosP) - WeylQntz(mSinQ, mSinP));

%s = -1i/sH.hbar * tau;
%matH = matT + matV;
%matH = matT + matV + s/2 * (matT*matV - matV*matT);
[pevecs, pevalsmat] = eig(Phi);
[pevals, sindex] = sort(real(diag(pevalsmat)), 'descend');
pevecs = pevecs(:, sindex);
pstates= eigs2states(sH, pevecs, pevals);
utils.saveeigs(sH, pevals, pstates, 'basis', 'q', 'header', 'ham', 'savedir', 'Suris');

toc 
disp("const matU")
tic
matU = sU.SImatrix(T, V, 'tau', tau, 'order', siorder);
toc

tic
disp("eig. matU")
[uevecs, uevals] = eig(matU);
toc

tic
disp("sort matU")
sindex = sU.sortbynorm(uevecs, pevecs);
uevecs = uevecs(:, sindex);
uevals = diag(uevals(sindex, sindex) );
toc

ustates = eigs2states(sU, uevecs, uevals);
utils.saveeigs(sU, uevals, ustates, 'basis', 'q', 'header', 'uni', 'savedir', 'Suris');

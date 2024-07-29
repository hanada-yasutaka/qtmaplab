clear all
%private_addpath('Advanpix');
addpath("~/Dropbox/Packages/qtmaplab/");
%addpath('/nfs/qtmaplab/') %% path to qtmaplab
%addpath('~/Applications/Advanpix/') %% path to Advanpix for multiple precision arthmetics

%mp.Test()
%return 
k = (15/16)^2;
T = @(x) x.^2/2;
V = @(x) k*cos(x);
dT = @(x) x;
dV = @(x) -k*sin(x);

%mp.Digits(32); % 任意精度の仮数部の桁数

domain = [-2*pi 2*pi;-pi pi];

basis = 'p';

savedir = sprintf("k%f", k);
if ~exist(savedir, 'dir')
    mkdir(savedir);
end
path = sprintf("%s/splitting_n1.dat", savedir);

if exist(path, 'file')
    file = fopen(path, 'a');
else
    file = fopen(path, 'w');
    fprintf(file, "#dim,\t1/hbar,\t Delta E,\t Delta E_0\n");
end


%for dim = 1200:200:5000
dim = 50;
    disp(dim);
    %sH = SplitHamiltonian(dim, domain, basis);
    [sU, state] = SplitUnitary(dim, domain, basis);
    hbar = sU.hbar;

    siorder = 2;
    tau = 1; % 任意精度計算の場合はmpで評価してください

    disp("construct matU")
    tic
    matU = sU.SImatrix(T, V, 'tau', tau, 'order', siorder);
    toc

    tic
    disp("eig. matU")
    [uevecs, uevals] = eig(matU);
    toc
    cs1 = state.coherent(-pi, 0);
    cs2 = state.coherent( pi, 0);
    dcs1 = (cs1 + cs2)/sqrt(2);
    dcs2 = (cs1 - cs2)/sqrt(2);

    norms= dcs1.inner(uevecs);
    [M, idx0] = max(abs2(norms));
    
    norms= dcs2.inner(uevecs);
    [M, idx1] = max(abs2(norms));

    
    
    eigphase = angle(uevals);
    quasiene  = -hbar*eigphase;
    splitting1 = abs(quasiene(idx0) -quasiene(idx1));
    disp(splitting1)
    fprintf(file, "%d\t%f\t%s\n", dim, 1/hbar, splitting1);

    ustates = eigs2states(sU, uevecs, uevals);
    s1 = state.tostate(uevecs(:, idx0));
    s2 = state.tostate(uevecs(:, idx1));
    utils.savestate(sU,s1.qrep(), 'uni_evecs_q1.dat', 'basis', 'q', 'savedir', savedir);    
    utils.savestate(sU,s2.qrep(), 'uni_evecs_q2.dat', 'basis', 'q', 'savedir', savedir);        
    
    %ustates = eigs2states(sU, uevecs, uevals);

%end

fclose(file);

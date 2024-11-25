clear all
%private_addpath('Advanpix');
addpath("../"); %% path to qtmaplab
%addpath('~/Applications/Advanpix/') %% path to Advanpix for multiple precision arthmetics


%%% system setting 

%dim = 120;
dim = 50;
%mp.Digits(200);
%domain = mp('[-4*pi 4*pi;-15 15]');
domain = [-4*pi 4*pi;-15 15];

basis = 'p';

T = @(x) x.^2/2;
V = @(x) x.^2/2 - 2*cos(x);
dT = @(x) x;
dV = @(x) x + 2*sin(x);

%tau = mp('0.5');
tau = 0.3;

%%% BCH Hamiltonian setting
sH = SplitHamiltonian(dim, domain, basis);
matT = sH.matT(T);
matV = sH.matV(V);
s = -1j/sH.hbar*tau;

%%%% 3rd order BCH Hamiltonian %%%
order = 'VT'; % 'TV' or 'VT';
bch = BCH(matT, matV, s, order);
refH = bch.Hamiltonian(3);

[refevecs, refevalsmat] = eig(refH);
[refevals, sindex] = sort(real(diag(refevalsmat)));

refevecs = refevecs(:,sindex);
states = eigs2states(sH, refevecs, refevals);

%%%% Higher order BCH Hamiltonian %%%%
matH = bch.Hamiltonian(7);
[hevecs, hevals] = eig(matH);
sindex = sortindex(hevecs, refevecs); %% sorted by 3rd order BCH states

hevecs = hevecs(:, sindex);
hevals = diag(hevals(sindex, sindex) );
hstates = eigs2states(sH, hevecs, hevecs);

%%%% quantum map %%%%
[sU, state] = SplitUnitary(dim, domain, basis);
siorder = -1;
%QSIevolve = sU.SIevolve(T, V, 'tau', tau, 'order', siorder);
matU = sU.SImatrix(T, V, 'tau', tau, 'order', siorder);
[uevecs, uevals] = eig(matU);
sindex = sortindex(uevecs, refevecs); %% sorted by 3rd order BCH states

uevecs = uevecs(:, sindex);
uevals = diag(uevals(sindex, sindex) );
ustates = eigs2states(sU, uevecs, uevals);



CSI = SimplecticIntegrator(dT, dV, 'tau', tau, 'order', siorder);
sample = 100;
tmax = 500;
twopi = 2*pi;
q0 = linspace(domain(1,1), domain(1,2) , sample);
p0 = zeros(1, sample);
q1 = zeros(1, sample);
p1 = linspace(domain(2,1), domain(2,2), sample);
x = [[q0, q1]; [p0, p1]];
traj = [[]; []];

for i=1:tmax
    x = CSI.evolve(x);
%    q = x(1,:);
%    p = x(2,:);
%    q = q - floor((q - pi)/twopi)*twopi - twopi;
    traj = horzcat(traj, x);        
end    

%s = state.coherent(0, 0);
%s2 = state.coherent(mp('+pi'),mp('0'));
%s = hstates(1);

fig = figure('Position', [10 10 700 700]);
ax1 = axes('Position',[0.1  0.55  .38 .38],'Box','on');
%ax2 = axes('Position',[0.55 0.55 .38 .38],'Box','on');
ax3 = axes('Position',[0.1  0.09  .38 .38],'Box','on');
ax4 = axes('Position',[0.55 0.09 .38 .38],'Box','on');
axs = [ax1 ax1 ax3 ax4];

for i=1:dim
    sh = hstates(i); 
    su = ustates(i);     
        
    for ax=axs
        hold(ax, 'on');
        grid(ax, 'on');
    end
    
    %%% plot axs(1): qrep
    ax = axs(1);
    plot(ax, sh.q, log10(abs2( sh.qrep() )), '-', 'LineWidth', 3)
    plot(ax, su.q, log10(abs2( su.qrep() )), '-', 'LineWidth', 3)    
    title(ax, sprintf("norm=%f, n=%d-th eig.", norm(su.qrep()), i ));
    xlabel(ax, '$q$', 'Interpreter', 'latex', 'FontSize', 15);
    ylabel(ax, '$|\langle q|\psi_n\rangle|^2$', 'Interpreter', 'latex', 'FontSize', 15);
    
    %%% plot axs(3): hsmrep
    ax = axs(3);
    [x,y,z] = su.hsmrep();
    contour(ax, x, y, z, 10 , 'LineColor', 'none', 'Fill','on');
    d = scatter(ax, traj(1,:), traj(2,:), 1, '.');
    
    zmax = max(z, [], 'all');
    colormap(ax, flipud(hot));
    caxis(ax, [-inf zmax]) % colorbar scale
    cb = colorbar(ax,'westoutside');
    cb.Position = cb.Position - [0.12, 0, 0, 0]; % position of colorbar
    cb.Ticks=[];  % remove colorbar ticks
    xlabel(ax, '$q$', 'Interpreter', 'latex', 'FontSize', 15);
    ylabel(ax, '$p$', 'Interpreter', 'latex', 'FontSize', 15);
    axis(ax, reshape(double(domain.'), 1, []));
    
    %axis(ax, [-20 20 -10 10]);
    %%% plot axs(4): prep
    ax = axs(4);
    plot(ax, log10(abs2( sh.prep() )), sh.p, '-','LineWidth', 3)
    plot(ax, log10(abs2( su.prep() )), su.p, '-','LineWidth', 3)    
    title(ax, sprintf("norm=%f", norm(su.prep()) ))
    xlabel(ax, '$|\langle p|\psi_n\rangle|^2$', 'Interpreter', 'latex', 'FontSize', 15);
    ylabel(ax, '$p$', 'Interpreter', 'latex', 'FontSize', 15);
    
    linkaxes([axs(1),axs(3)], 'x');
    linkaxes([axs(3),axs(4)], 'y');
    %fprintf("%d\n", i);
    
    savedir="test";    
    if ~exist(savedir, 'dir')
       mkdir(savedir);
    end
    saveas(gcf, sprintf("test/test_%d.png", i));
    
    for ax=axs
        hold(ax, 'off')
    end
    
    fprintf("press button to the next:\n");
    waitforbuttonpress
        
    colorbar(axs(3), 'off');
    for ax=axs
        cla(ax);
    end
    
end


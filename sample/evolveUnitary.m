clear all
private_addpath('Advanpix/');

dim = 100;
%mp.Digits(150);
domain = mp('[-pi pi;-pi pi]');
%domain = [0 pi;0 pi];

k = 1;
tau = 1; %mp('0.1');
%tau = 0.1
basis = 'p';

T = @(x) x.^2/2;
V = @(x) cos(x);
dT = @(x) x;
dV = @(x) -sin(x);

sH = SplitHamiltonian(dim, domain, basis);
[sU, state] = SplitUnitary(dim, domain, basis);
siorder = 2;
QSIevolve = sU.SIevolve(T, V, 'tau', tau, 'order', siorder);
CSI = SimplecticIntegrator(dT, dV, 'tau', tau, 'order', siorder)


matT = sH.matT(T);
matV = sH.matV(V);
matH = matT + matV;
%[hevecs, hevalsmat] = eig(matH);
%[hevals, sindex] = sort(real(diag(hevalsmat)));
%hevecs = hevecs(:, sindex);
%hstates= eigs2states(sH, hevecs, hevals);

sample = 20;
tmax = 300;
twopi = 2*pi;
q0 = linspace(domain(1,1), domain(1,2) , sample);
p0 = zeros(1, sample);
q1 = zeros(1, sample);
p1 = linspace(domain(2,1), domain(2,2), sample);
x = [[q0, q1]; [p0, p1]];
traj = [[]; []];

for i=1:tmax
    x = CSI.evolve(x);
    q = x(1,:);
    p = x(2,:);
    q = q - floor((q - pi)/twopi)*twopi - twopi;
    traj = horzcat(traj, [q;p]);        
end    

s = state.coherent(-pi, 0);
%s2 = state.coherent(mp('+pi'),mp('0'));
%s = hstates(1);

fig = figure('Position', [10 10 700 700]);
ax1 = axes('Position',[0.1  0.55  .38 .38],'Box','on');
%ax2 = axes('Position',[0.55 0.55 .38 .38],'Box','on');
ax3 = axes('Position',[0.1  0.09  .38 .38],'Box','on');
ax4 = axes('Position',[0.55 0.09 .38 .38],'Box','on');
axs = [ax1 ax1 ax3 ax4];

for i=0:100
        
    for ax=axs
        hold(ax, 'on');
    end
    
    %%% plot axs(1): qrep
    ax = axs(1);
    plot(ax, s.q, log10(abs2( s.qrep() )), '-')
    title(ax, sprintf("norm=%f", norm(s.qrep()) ));
    xlabel(ax, '$q$', 'Interpreter', 'latex', 'FontSize', 15);
    ylabel(ax, '$|\langle q|\psi_n\rangle|^2$', 'Interpreter', 'latex', 'FontSize', 15);
    
    %%% plot axs(3): hsmrep
    ax = axs(3);
    [x,y,z] = s.hsmrep('ismp', true);
    contour(ax, x, y, log10(z), 10 , 'LineColor', 'none', 'Fill','on');
    d = scatter(ax, traj(1,:), traj(2,:), 1, '.');
    
    zmax = max(z, [], 'all');
    colormap(ax, flipud(hot));
    caxis(ax, double([-inf zmax])) % colorbar scale
    cb = colorbar(ax,'westoutside');
    cb.Position = cb.Position - [0.12, 0, 0, 0]; % position of colorbar
    cb.Ticks=[];  % remove colorbar ticks
    xlabel(ax, '$q$', 'Interpreter', 'latex', 'FontSize', 15);
    ylabel(ax, '$p$', 'Interpreter', 'latex', 'FontSize', 15);
    axis(ax, reshape(double(domain).', 1, []));
    
    %%% plot axs(4): prep
    ax = axs(4);
    plot(ax, log10(abs2( s.prep() )), s.p, '-')
    title(ax, sprintf("norm=%f", norm(s.prep()) ))
    xlabel(ax, '$|\langle p|\psi_n\rangle|^2$', 'Interpreter', 'latex', 'FontSize', 15);
    ylabel(ax, '$p$', 'Interpreter', 'latex', 'FontSize', 15);
    
    %fprintf("%d\n", i);
    %savefig(sprintf("test/test_%d.png", i));
    for ax=axs
        hold(ax, 'off')
    end
    
    fprintf("press button to the next:\n");
    waitforbuttonpress
        
    colorbar(axs(3), 'off');
    for ax=axs
        cla(ax);
    end
    
    utils.savestate(sH, s, sprintf('evolve_t%d.dat', i), 'basis', 'q', 'savedir', 'Data');    
    s = QSIevolve(s);
end


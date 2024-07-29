clear all
%private_addpath('Advanpix');
addpath("~/Dropbox/Packages/qtmaplab/"); %% path to qtmaplab
%addpath('~/Applications/Advanpix/') %% path to Advanpix for multiple precision arthmetics

%mp.Digits(100); % 仮数部の桁数
dim = 50; 

%domain = mp('[-pi pi;-pi pi]');
domain = [0 2*pi; -1/2*pi 1/2*pi]

basis = 'p';
sH = SplitHamiltonian(dim, domain, basis);
sU = SplitUnitary(dim, domain, basis);

T = @(x) x.^2/2;
V = @(x) cos(x);
dT = @(x) x;
dV = @(x) -sin(x);

matT = sH.matT(T);
matV = sH.matV(V);

siorder = 1;
%tau = mp('1');
tau = 1;

tic 
disp("eig Hamiltonian")
s = -1i/sH.hbar * tau;
%matH = matT + matV;
matH = matT + matV + s/2 * (matT*matV - matV*matT);
[hevecs, hevalsmat] = eig(matH);
[hevals, sindex] = sort(real(diag(hevalsmat)));
hevecs = hevecs(:, sindex);
hstates= eigs2states(sH, hevecs, hevals);
%saveobj(hstates(1), 'abc.mat')
utils.saveeigs(sH, hevals, hstates, 'basis', 'q', 'header', 'ham', 'savedir', 'Data');

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
sindex = sU.sortbynorm(uevecs, hevecs);
uevecs = uevecs(:, sindex);
uevals = diag(uevals(sindex, sindex) );
toc

ustates = eigs2states(sU, uevecs, uevals);
utils.saveeigs(sU, uevals, ustates, 'basis', 'q', 'header', 'uni', 'savedir', 'Data');


CSI = SimplecticIntegrator(dT, dV, 'tau', double(tau), 'order', siorder);

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
    %q = q - floor((q - pi)/twopi)*twopi - twopi;
    traj = horzcat(traj, [q;p]);        
end    

fig = figure('Position', [10 10 700 700]);
ax1 = axes('Position',[0.1  0.55  .38 .38],'Box','on');
ax2 = axes('Position',[0.55 0.55 .38 .38],'Box','on');
ax3 = axes('Position',[0.1  0.09  .38 .38],'Box','on');
ax4 = axes('Position',[0.55 0.09 .38 .38],'Box','on');
axs = [ax1 ax2 ax3 ax4];


for i=1:dim
    s = ustates(i); 

    for ax=axs
        hold(ax, 'on');
    end
        
    %%% plot axs(1): qrep
    ax = axs(1);
    plot(ax, s.q, log10(abs2( s.qrep() )), '-')
    title(ax, sprintf("norm=%f", norm(s.qrep()) ));
    xlabel(ax, '$q$', 'Interpreter', 'latex', 'FontSize', 15);    
    ylabel(ax, '$|\langle q|\psi_n\rangle|^2$', 'Interpreter', 'latex', 'FontSize', 15);
    
    
    %%% plot axs(2): eigenvalues 
    ax = axs(2);
    theta = linspace(-pi, pi, 100);
    z = exp(1j*theta);
    plot(ax, real(z), imag(z), '-k');
    scatter(ax, real(uevals), imag(uevals) );
    scatter(ax, real(s.eigenvalue), imag(s.eigenvalue),100, 'filled');
    title(ax, sprintf("%d-th eigs, $u_n$=%f%+fi", i, real(s.eigenvalue), imag(s.eigenvalue)), 'Interpreter', 'latex');    
    xlabel(ax, '$\mathrm{Re}(u_n)$', 'Interpreter', 'latex', 'FontSize', 15);
    ylabel(ax, '$\mathrm{Im}(u_m)$', 'Interpreter', 'latex', 'FontSize', 15);    
    
    
    %%% plot axs(3): hsmplot
    ax = axs(3);
    ds = double(s);
    [x,y,z] = ds.hsmrep('periodic', true);
    contour(ax, x, y, z, 10 , 'LineColor', 'none', 'Fill','on');
    
    d = scatter(ax, traj(1,:), traj(2,:), 1, '.');%, 'MarkerSize', 1, 'Marker', 'o')    
    zmax = max(z, [], 'all');
    colormap(ax, flipud(hot));
    %caxis(ax, [-inf zmax]) % colorbar scale
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
    
    for ax=axs
        colorbar('off');
        cla(ax);        
    end    
end





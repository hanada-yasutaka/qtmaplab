clear all
private_addpath('Advanpix');

dim = 10;
domain = mp('[-pi pi;-2*pi 2*pi]');
%domain = [-pi pi;-pi pi];
basis = 'p';
sH = SplitHamiltonian(dim, domain, basis);

T = @(x) x.^2/2;
V = @(x) cos(x);
matT = sH.matT(T);
matV = sH.matV(V);

matH = matT + matV;

[evecs, evalsmat] = eig(matH);
[evals, sindex] = sort(real(diag(evalsmat)));

evecs = evecs(:,sindex);
states = eigs2states(sH, evecs, evals);
utils.saveeigs(sH, evals, states, 'basis', 'q', 'header', 'ham', 'savedir', 'Data');
%save 'Data/matH.mat' matH;
%save 'Data/states.mat' states;

sample = 100;
q = linspace(domain(1,1), domain(1,2), sample);
p = linspace(domain(2,1), domain(2,2), sample);
[Q, P] = meshgrid(q, p);
H = T(P) + V(Q);


fig = figure('Position', [10 10 700 700]);
ax1 = axes('Position',[0.1  0.55  .38 .38],'Box','on');
ax2 = axes('Position',[0.55 0.55 .38 .38],'Box','on');
ax3 = axes('Position',[0.1  0.09  .38 .38],'Box','on');
ax4 = axes('Position',[0.55 0.09 .38 .38],'Box','on');
axs = [ax1 ax2 ax3 ax4];

for i=1:dim
    s = states(i);     

    for ax=axs
        hold(ax, 'on');
    end
   
    %%% plot axs(1); qrep
    ax = axs(1);
    plot(ax, s.q, log10( abs2(s.qrep() )) );
    title(ax, sprintf("norm=%f", norm(s.qrep())) );    
    xlabel(ax, '$q$', 'Interpreter', 'latex', 'FontSize', 15);    
    ylabel(ax, '$|\langle q|\psi_n\rangle|^2$', 'Interpreter', 'latex', 'FontSize', 15);
    
    %%% plot axs(2): eigen energy
    ax = axs(2);
    scatter(ax, 1:dim, evals);
    scatter(ax, i, s.eigenvalue, 100,'filled')
    title(ax, sprintf("%d-th eigs,E_n=%f", i, s.eigenvalue));
    xlabel(ax, '$n$', 'Interpreter', 'latex', 'FontSize', 15);
    ylabel(ax, '$E_n$', 'Interpreter', 'latex', 'FontSize', 15);    

    %%% plot axs(3): hsmplot
    ax = axs(3);
    [x,y,z] = s.hsmrep('gridnum', 50, 'periodic', true);    
    contour(ax, Q, P, H, 10, 'LineColor', [0 0 0],'LineWidth', 0.5);        
    contour(ax, x, y, z, 10 , 'LineColor', 'none', 'Fill','on');
    contour(ax, Q, P, H, [1 1]*s.eigenvalue, 'LineColor', [0 1 1],'LineWidth', 2);                
    class(H)
    [1 1]*s.eigenvalue
    zmax = max(z, [], 'all');
    colormap(ax, flipud(hot));
    caxis(ax, [-inf double(zmax)]) % colorbar scale
    cb = colorbar(ax,'westoutside');
    cb.Position = cb.Position - [0.12, 0, 0, 0]; % position of colorbar
    cb.Ticks=[];  % remove colorbar ticks
    xlabel(ax, '$q$', 'Interpreter', 'latex', 'FontSize', 15);
    ylabel(ax, '$p$', 'Interpreter', 'latex', 'FontSize', 15);    
    
    %%% plot axs(4): prep
    ax = axs(4);
    plot(ax, log10( abs2(s.prep()) ), s.p)
    title(ax, sprintf("norm=%f", norm(s.prep())) );
    xlabel(ax, '$|\langle p|\psi_n\rangle|^2$', 'Interpreter', 'latex', 'FontSize', 15);
    ylabel(ax, '$p$', 'Interpreter', 'latex', 'FontSize', 15);        
    
    
    for ax=axs
        hold(ax, 'off')
    end
    
    fprintf("press/click to the next:\n");
    waitforbuttonpress
    
    for ax=axs
        colorbar('off');
        cla(ax);        
    end
    
end

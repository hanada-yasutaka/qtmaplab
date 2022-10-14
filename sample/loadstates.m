private_addpath('Advanpix');

list = dir("Data/ham_qrep*.dat");
states = utils.loadeigs(list);

fig = figure('Position', [10 10 700 700]);
ax1 = axes('Position',[0.1  0.55  .38 .38],'Box','on');
%ax2 = axes('Position',[0.55 0.55 .38 .38],'Box','on');
ax3 = axes('Position',[0.1  0.09  .38 .38],'Box','on');
ax4 = axes('Position',[0.55 0.09 .38 .38],'Box','on');
axs = [ax1 ax1 ax3 ax4];

for i = 1:length(states)

    for ax=axs
        hold(ax, 'on');
        grid(ax, 'on');
    end
    
    s = states(i);
    %%% plot axs(1): qrep
    ax = axs(1);
    plot(ax, s.q, log10(abs2( s.qrep() )), '-', 'LineWidth', 3)
    title(ax, sprintf("norm=%f, n=%d", norm(s.qrep()), i ));
    xlabel(ax, '$q$', 'Interpreter', 'latex', 'FontSize', 15);
    ylabel(ax, '$|\langle q|\psi_n\rangle|^2$', 'Interpreter', 'latex', 'FontSize', 15);
    
    %%% plot axs(3): hsmrep
    ax = axs(3);
    [x,y,z] = s.hsmrep();
    contour(ax, x, y, log10(z), 10 , 'LineColor', 'none', 'Fill','on');
    %d = scatter(ax, traj(1,:), traj(2,:), 1, '.');
    
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
    plot(ax, log10(abs2( s.prep() )), s.p, '-','LineWidth', 3)
    title(ax, sprintf("norm=%f", norm(s.prep()) ))
    xlabel(ax, '$|\langle p|\psi_n\rangle|^2$', 'Interpreter', 'latex', 'FontSize', 15);
    ylabel(ax, '$p$', 'Interpreter', 'latex', 'FontSize', 15);
    
    linkaxes([axs(1),axs(3)], 'x');
    linkaxes([axs(3),axs(4)], 'y');
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
 
end

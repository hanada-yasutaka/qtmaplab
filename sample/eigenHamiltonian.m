clear all
AdvanpixMCT='AdvanpixMCT-4.8.3.14460/';
if ismac
    addpath("/Users/hanada/OneDrive/Packages/qtmaplab/");
    addpath(sprintf("/Users/hanada/Applications/%s",AdvanpixMCT));
elseif isunix
    [~, name] = system('hostname');
    if strcmp(strtrim(name), 'bohigas')
        addpath("/home/hanada/OneDrive/Packages/qtmaplab/");
        addpath(sprintf("/home/hanada/Applications/%s",AdvanpixMCT));        
    else
        addpath("/nfs/qtmaplab/");
        addpath("/nfs/AdvanpixMCT-4.8.3.14440/");
    end
end


dim = 50;
%domain = mp('[-pi pi;-2*pi 2*pi]');
%domain = [-pi pi;0 2*pi];
%domain = [0 2*pi;0 2*pi];
%domain = [-pi pi;-pi pi];
domain = [0 2*pi;-pi pi];

%domain = [-2*pi 2*pi;-2*pi 2*pi];
basis = 'p';
sH = SplitHamiltonian(dim, domain, basis);

T = @(x) x.^2/2;
V = @(x) cos(x);
matT = sH.matT(T);
matV = sH.matV(V);
%class(T(1))
%assertWarningFree(@T(1))
%assertWarningFree(T(1))
%matT = sH.matT(T);
%matV = sH.matV(V);
%return 
matH = matT + matV;
[evecs, evalsmat] = eig(matH);
[evals, sindex] = sort(real(diag(evalsmat)));

evecs = evecs(:,sindex);
states = eigs2states(sH, evecs, evals);

sample = 100;
q = linspace(domain(1,1), domain(1,2), sample);
p = linspace(domain(2,1), domain(2,2), sample);
[Q, P] = meshgrid(q, p);
H = T(P) + V(Q);

fig = figure();
axs = [subplot(2,2,1) subplot(2,2,2) subplot(2,2,3) subplot(2,2,4) ] ;


for i=1:dim
    s = states(i);     
    
    for ax=axs
        hold(ax, 'on');
    end
   

    plot(axs(1), s.q, log10( abs2(s.qrep() )) );
    title( sprintf("norm=%f", norm(s.qrep())) );    
    
    scatter(axs(2), 0:dim-1, evals);
    scatter(axs(2), i, s.eigenvalue, 100,'filled')
    title(sprintf("%d-th eigs,E_n=%f", i, s.eigenvalue));

    
    contour(axs(3), Q, P, H, 10, 'LineColor', [0 0 0],'LineWidth', 0.5);
    [x,y,z] = s.hsmrep(50);    
    contour(axs(3), x, y, z, 20, 'LineColor', 'flat','LineWidth', 2);
    zmax = max(z, [], 'all');
    
    caxis([0 zmax])
    cb = colorbar(axs(3),'westoutside')
    cb.Position = cb.Position - [0.12, 0, 0, 0]

    plot(axs(4), log10( abs2(s.prep()) ), s.p)
    title(sprintf("norm=%f", norm(s.prep()) ));
    for ax=axs
        hold(ax, 'off')
    end
    fprintf("press to the next:\n");
    waitforbuttonpress
    for ax=axs
        cla(ax);
    end
end

function y = fT(x)
y = x.^2/2;
end

function y = fV(x)
y = cos(x);
end
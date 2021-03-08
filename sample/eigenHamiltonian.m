clear all
%addpath("/nfs/AdvanpixMCT-4.8.3.14440/");
%addpath("/Users/hanada/Documents/MATLAB/qtmaplab/");
addpath("/Users/hanada/OneDrive/Packages/qtmaplab/")

dim = 50;
domain = [-pi pi;0 2*pi];
basis = 'p';
sH = SplitHamiltonian(dim, domain, basis);

%V = @(x) x .^2 / 2;
V = @(x) sin(x+2);
T = @(x) sin(x+1);

matT = sH.matT(T);
matV = sH.matV(V);

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

for i=1:dim
    s = states(i);     
    tiledlayout(2,2)
    nexttile;
    str = sprintf("norm=%f", norm( s.qrep() ) );
    plot(s.q, log10( abs2(s.qrep() )) );
    title(str);    
    
    nexttile;
    e = s.eigenvalue;
    ns = linspace(0,dim-1,dim);
    scatter(ns, evals);
    hold on
    scatter(i, s.eigenvalue,100,'filled')
    str = sprintf("%d-th eigs,E_n=%f", i, e); 
    title(str);
    
    hold off
    
    
    nexttile;
    
    contour(Q, P, H, 10, 'LineColor', [0 0 0],'LineWidth', 0.5);
    hold on
    s = states(i); 
    [x,y,z] = s.hsmrep();
    contour(x, y, z, 20, 'LineColor', 'flat','LineWidth', 2);
    zmax = max(z, [], 'all');
    caxis([0 zmax])
    colorbar
    hold off
    
    nexttile;
    plot(log10( abs2(s.prep()) ), s.p)
    title(sprintf("norm=%f", norm(s.prep()) ));
    

    fprintf("press enter to the next:\n")
    waitforbuttonpress
end
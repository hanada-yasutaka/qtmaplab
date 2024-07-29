clear all
%private_addpath('Advanpix');
addpath("~/Dropbox/Packages/qtmaplab/"); %% path to qtmaplab
%addpath('~/Applications/Advanpix/') %% path to Advanpix for multiple precision arthmetics

domain = mp('[-2*pi 2*pi;-2*pi 2*pi]');
%domain = [-2*pi 2*pi;-2*pi 2*pi];
basis = 'p';

T = @(x) x.^2/2;
V = @(x) cos(x);

of = fopen("test_spl_mp_apple.dat", "w");

sample = 100;
q = linspace(domain(1,1), domain(1,2), sample);
p = linspace(domain(2,1), domain(2,2), sample);
[Q, P] = meshgrid(q, p);
H = T(P) + V(Q);

tic 
for dim = 2:1:200
    sH = SplitHamiltonian(dim, domain, basis);
    matT = sH.matT(T);
    matV = sH.matV(V);

    matH = matT + matV;
    [evecs, evalsmat] = eig(matH);

    [evals, sindex] = sort(real(diag(evalsmat)));
    evecs = evecs(:,sindex);
    
    if mod(dim,2) == 0    
        evals0 = evals;
    else
        evals1 = evals;
    end
    spl = evals(2) - evals(1);
    evec0 = evecs(:,sindex);
    fprintf("dim=%d\n", dim);
    fprintf(of, "%d %.18e\n", dim, spl);
    states = eigs2states(sH, evecs, evals);

%     s = states(2);
%     tiledlayout(2,2)
%     nexttile;
%     str = sprintf("dim=%d norm=%f", dim,norm( s.qrep() ) );
%     plot(s.q, log10( abs2(s.qrep()) ) );
%     ylim([-30 1]);
%     title(str);    
%     
%     nexttile;    
%     if mod(dim,2)==1
%         scatter(1:length(evals0), evals0);
%         hold on
%         scatter(1:length(evals1), evals1);
%         xlim([0 10]);
%         ylim([-1 3]); 
%         hold off
%         fprintf("%.18e\n%.18e\n", evals0(2) - evals0(1), evals(2) - evals(1));
%     end
%     
%     nexttile;
%     
%     contour(Q, P, H, 10, 'LineColor', [0 0 0],'LineWidth', 0.5);
%     hold on
%     [x,y,z] = s.hsmrep(100);
%     contour(x, y, z, 20, 'LineColor', 'flat','LineWidth', 2);
%     zmax = max(z, [], 'all');
%     caxis([0 zmax])
%     colorbar
%     contour(Q, P, H, [s.eigenvalue, 1], 'LineColor', 'r','LineWidth', 3);    
%     
%     hold off
%     
%     nexttile;
%     plot(log10( abs2(s.prep()) ), s.p)
%     xlim([-30 1]);
%     title(sprintf("norm=%f", norm(s.prep()) ));
%     
% 
%     fprintf("press enter to the next:\n")
%     waitforbuttonpress
end
toc
fclose(of);
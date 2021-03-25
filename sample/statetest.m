clear all
addpath("/Users/hanada/OneDrive/Packages/qtmaplab/")

dim = 100;
%domain = [-pi pi; -pi pi];
domain = [0 2*pi; 0 2*pi];
basis = 'p';
sysinfo=SystemInfo(dim, domain);
state = FundamentalState(sysinfo, basis);

s = state.coherent(2, 2);
if strcmp(basis, 'q')
    ss = s.q2p();
    %ss = ss.p2q()
    %ss = ss.q2p()
else
    ss = s.p2q();
    %ss = ss.p2q()
end

%vec = state.delta(-2);
%vec = vec.p2q();

tiledlayout(2,2)
nexttile;
plot(s.q, log10(abs2( s.qrep() )), '-', 'LineWidth', 3)
xlim(domain(1,:))
hold on 
plot(ss.q, log10(abs2( ss.qrep() )), '-', 'LineWidth', 2)
hold off
title(sprintf("norm=%f,norm=%f", norm(s.qrep())), norm(ss.qrep() ) )
nexttile;
[x,y,z] = s.hsmrep();
contour(x, y, z, 10, 'LineWidth', 3);
axis(reshape(domain.', 1, []));

nexttile;

[x,y,z] = ss.hsmrep();
contour(x, y, z, 10, 'LineWidth', 3);
axis(reshape(domain.', 1, []));

nexttile
plot(log10(abs2( s.prep() )), s.p, '-', 'LineWidth', 3)
ylim(domain(2,:))
hold on 
plot(log10(abs2( ss.prep() )), ss.p, '-', 'LineWidth', 2)
hold off
title(sprintf("norm=%f,norm=%f", norm(s.prep())), norm(ss.prep() ) )


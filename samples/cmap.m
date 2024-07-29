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
        addpath(sprintf("/nfs/%s", AdvanpixMCT));
    end
end


T = @(x) x.^2/2;
V = @(x) cos(x);

dT = @(x) x;
dV = @(x) -sin(x);

CSI = SimplecticIntegrator(dT, dV, 1, 1);
sample = 10;
q0 = linspace(0, 2*pi , sample);
p0 = zeros(1, sample);
q1 = zeros(1, sample);
p1 = linspace(-2*pi, 2*pi, sample);

x = [[q0, q1]; [p0,p1]];
traj = [[]; []];


tmax = 1000;
twopi = 2*pi;
%traj = CSI.getTraj(x, tmax);
for i=1:tmax
    x = CSI.evolve(x);
    traj = horzcat(traj, x);
end    

q = traj(1,:);
p = traj(2,:);
q = q - floor((q - pi)/twopi)*twopi - twopi;
%scatter(traj(1,:), traj(2,:), '.')
scatter(q, p, '.')
return 




%[q, p] = [q - obj.dV(q) * tau * c; p];
%[q, p] = [q + obj.dT(p) * tau * c; p];

%x = CSI.evolveV(x);
q = x(1, :);
p = x(2, :);

%xx = [xx(1,:); xx(1,:) + dT(xx(2,:))]
p = p - dV(q);
q = q;
q = q + dT(p);

plot(q, p, '-')
hold on 


%qq = xx(1, :);
%pp = xx(2, :);
%xx = [pp - dV(qq); qq];
xx = [x(1,:); x(2,:)];
xx = CSI.evolveT(xx);
xx = CSI.evolveV(xx);
%xx = [q0; p0];
plot(xx(1,:), xx(2,:), 'x')

xx = [x(1,:); x(2,:)];
q = xx(1,:);
p = xx(2,:);
xx = [q; p];
xx = CSI.evolveTV(xx)
%[q,p] = [q; p - dV(q)]
%[q,p] = [q + dT(p); p]

%scatter(q, p,"o")
%plot(q0, p0 - dV(q0))
%scatter(x(1,:), x(2,:))

hold off
%plot(q, p)
%for i=1:tmax

    %q = q - floor((q - pi)/twopi)*twopi - twopi;
    %p0 = p0 - floor((p0 - pi)/twopi)*twopi - twopi;
    %q0 = q0 - floor(q0/twopi)*twopi;
    %p0 = p0 - floor(p0/twopi)*twopi;    
%    trajq = [trajq, x(1)];
%    trajp = [trajp, x(2)];
%end    

%scatter(trajq, trajp, '.')

%axs = [subplot(2,2,1) subplot(2,2,2) subplot(2,2,3) subplot(2,2,4) ] ;


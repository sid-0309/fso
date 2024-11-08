import fso.outdoor.ao

% plt_resolution = 101;
% dfs = linspace(0,0.5,plt_resolution);
% etas = zeros(plt_resolution,1);
% 
% mu = [0 0];
% sigma = [3e-6 0;0 3e-6];
% fn = normal_phase(x,y,mu,sigma);
% 
% for a = 1:plt_resolution
% 
%     sys = ao;
%     f = (12e-3)/dfs(a);
%     w0 = 8e-6;
%     w = (1550e-9)*(f)/(pi*w0);
% 
%     [x,y] = sys.receiver(sys.D/2 ,sys.D/2 ,101 ,101);
% 
%     phi = sys.received_wf(101,101,1,false);
% 
%     ea = sys.receiver_field(phi);
%     fa = sys.of_mode_field(x,y,w);
% 
%     eta = sys.coupling_efficiency(ea,fa);
%     etas(a) = eta*100;
% end
% 
% plot(dfs,etas);


%% Single run
%import fso.outdoor.ao
% sys = ao;
% df = input("Enter D/f: ");
% f = (12e-3)/df;
% w0 = 8e-6;
% w = (1550e-9)*(f)/(pi*w0);
% 
% mu = [0 0];
% sigma = [3e-6 0;0 3e-6];
% fn = normal_phase(x,y,mu,sigma);
% 
% [x,y] = sys.receiver(sys.D/2 ,sys.D/2 ,101 ,101);
% 
% phi = sys.received_wf(101,101,1,false);
% 
% ea = sys.receiver_field(phi);
% fa = sys.of_mode_field(x,y,w);
% 
% eta = sys.coupling_efficiency(ea,fa);
% 
% disp(eta*100);


%% Zernike Fitting

import fso.outdoor.ao
N = input("Enter N: ");
sys = ao;
[x,y] = sys.receiver(sys.D/2 ,sys.D/2 ,101 ,101);
mu = [0 0];
sigma = [3e-6 0;0 3e-6];
fn = normal_phase(x,y,mu,sigma);
fn_norm = fn./max(max(fn));
phi = sys.received_wf(101,101,fn_norm);

[coefficients, zernikeModes] = sys.fitZernikeModes(phi,N);

% Initialize the fitted Zernike surface
zernike_surface = zeros(size(phi));

% Loop through each Zernike mode and apply the coefficients
for k = 1:N
    zernike_surface = zernike_surface + coefficients(k) * zernikeModes{k};
end

%% Zernike Fitting Error Plot

import fso.outdoor.ao
% d = input("Receiver divisions : ");
d = 201;
% N = input("Enter max number of Zernike modes : ");
N = 40;
sys = ao;

% Define receiver
[x,y] = sys.receiver(sys.D/2,sys.D/2,d,d);

% For modeling received wavefront as gaussian
mu = [0 0];
sigma = [3e-6 0;0 3e-6];

% Log-Normal Noise
% mulog = input("mu : ");
% logNoise = lognrnd(0,mulog,[d,d]);
% fn = normal_phase(x,y,mu,sigma).*logNoise;

% Gamma-Gamma Noise
alpha = 1e3;
beta = 1e3;
g1 = gamrnd(alpha,1/alpha,[d,d]);
g2 = gamrnd(beta,1/beta,[d,d]);
fn = normal_phase(x,y,mu,sigma).*(g1.*g2);

% Adding AWGN
fn = fn./max(max(fn))+rand(d,d)*3e-2;
%Normalizing wavefront
fn_norm = fn./(max(max(fn)));

% figure();
% noise = g1.*g2;
% noise = noise./max(max(noise));
% surf(fn./max(max(fn))-noise);
% figure();
% surf(fn_norm);

% Vector for storing mean error, mean of abosute error and standard
% deviation respectively
m = zeros(1,N);
am = zeros(1,N);
s = zeros(1,N);

for i = 1:N
    [coefficients, zernikeModes] = sys.fitZernikeModes(fn_norm,i);
    zernike_surface = zeros(size(fn_norm));

    % finding the corrected wavefront
    for k = 1:i
        zernike_surface = zernike_surface + coefficients(k) * zernikeModes{k};
    end

    % flattening the wavefront for calculation
    err = reshape(zernike_surface-fn_norm,1,[]);
    am(i) = mean(abs(err));
    m(i) = mean(err);
    s(i) = std(err);
end
figure();
subplot(1,3,1);
% plot(1:N,m);
errorbar(1:N,m,-s,s);
subplot(1,3,2)
plot(1:N,s);
subplot(1,3,3);
plot(1:N,am);

%% Plot Zernike Surfaces

import fso.outdoor.ao
sys = ao;
rows = 201;
cols = 201;
[Zx, Zy] = meshgrid(linspace(-1, 1, cols), linspace(-1, 1, rows));
r = sqrt(Zx.^2 + Zy.^2);
mask = r <= 1;
theta = atan2(Zy, Zx);


n = input("Enter max radial mode : ");

surfs = cell(1,(n+1)*(n+2)/2);

pltn = 1;
k = 1;
for i = 0:n
    for j = -i:2:i
        surfs{k} = sys.zernike(r,theta,i,j).*mask;
        subplot(n+1,n+1,pltn);
        surf(surfs{k});
        shading interp
        k = k+1;
        pltn = pltn+1;

    end
    if mod(pltn,n+1) ~= 0
        pltn = pltn+(n+1-mod(pltn,n+1)+1);
    else
        pltn = pltn+1;
    end
end
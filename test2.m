import fso.outdoor.ao

% plt_resolution = 101;
% dfs = linspace(0,0.5,plt_resolution);
% etas = zeros(plt_resolution,1);
% 
% % mu = [0 0];
% % sigma = [3e-6 0;0 3e-6];
% % fn = normal_phase(x,y,mu,sigma,1e-4);
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
% sys = ao;
% df = input("Enter D/f: ");
% f = (12e-3)/df;
% w0 = 8e-6;
% w = (1550e-9)*(f)/(pi*w0);
% 
% mu = [0 0];
% sigma = [3e-6 0;0 3e-6];
% % fn = normal_phase(x,y,mu,sigma);
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
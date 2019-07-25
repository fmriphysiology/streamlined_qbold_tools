function [r2p, dbv, oef, dhb] = ase_qbold_3d(nii_name,start_tau,delta_tau,end_tau,nii_out)
% nii_name: Z-averaged ASE dataset
% start_tau [ms] value of tau for first volume
% delta_tau [ms] tau step size
% set end_tau [ms] ... value of tau for last volume
% nii_out ... output nifti's and save figs? 1=Yes,0=No 

%% Load dataset
[V,dims,scales,bpp,endian] = read_avw(nii_name);
[x y z v] = size(V);
%disp(size(V))

%% Fit R2'
% Constants
Hct = 0.34; % hct ratio in small vessels
dChi0 = 0.264*10^-6; % ppm, sus difference between fully oxy & deoxy rbc's
gamma = 2.675*10^4; % rads/(secs.Gauss) gyromagnetic ratio
B0 = 3*10^4; %Gauss, Field strength
phi = 1.34; % mlO2/gHb
k = 0.03; % conversion factor Hct (% rbc's in blood) to Hb (iron-containing molecules in the rbc's used to transport o2)


% X
% Convert tau to seconds
tau = [start_tau:delta_tau:end_tau];
tau = tau.*10^-3;

% Check that tau matches the number of volumes 

if (length(tau) ~= v)
    disp('List of Tau values doesn''t match the number of volumes') 
    %disp(tau)
    sprintf('Number of volumes = %1.0f', v)
else

% Y
ln_Sase = log(V); 
ln_Sase(isnan(ln_Sase)) = 0; 
ln_Sase(isinf(ln_Sase)) = 0;
ln_Saser = reshape(ln_Sase,x*y*z,v);

Tc = 0.015; % cutoff time for monoexponential regime [s]
tau_lineID = find(tau > Tc); % tau's to be used for R2' fitting  
s0_id = find(tau == 0); % find spin echo index
%w = 1./tau(tau_lineID)'; % weightings for lscov
p = zeros(x,y,z,2); 

%% LSCOV: fit linear regime
X = [ones(size(tau(tau_lineID)')) -tau(tau_lineID)' ones(size(tau(tau_lineID)'))];
X = [0 0 1; X];

Y = ln_Saser(:,tau_lineID)';
Y = [ln_Saser(:,s0_id)'; Y];

p = lscov(X,Y); %no weighting

range_r2p = [0 7];
range_dbv = [0 0.12];
range_oef = [0 1.0];
range_dhb = [0 10];

% Calculate Physiological Parameters
dbv = reshape(p(1,:),x,y,z);
r2p = reshape(p(2,:),x,y,z); 
s0 = reshape(p(3,:),x,y,z); 
oef = r2p./(dbv.*gamma.*(4./3).*pi.*dChi0.*Hct.*B0);
dhb = r2p./(dbv.*gamma.*(4./3).*pi.*dChi0.*B0.*k);

% Display parameter maps
imgH_r2p = brain_montage(r2p, range_r2p, 0, 'R_2''','R2''','JetBlack', '[s^{-1}]',[1 z]);
if nii_out, print('r2p.eps','-depsc2','-r300'), end
imgH_dbv = brain_montage(dbv, range_dbv, 0, 'DBV','DBV','JetBlack', '[%]',[1 z]);
if nii_out, print('dbv.eps','-depsc2','-r300'), end
imgH_oef = brain_montage(oef, range_oef, 0, 'OEF','OEF','JetBlack', '[%]',[1 z]);
if nii_out, print('oef.eps','-depsc2','-r300'), end
imgH_dhb = brain_montage(dhb, range_dhb, 0, 'dHb','dHb','JetBlack', '[g.dl^{-1}]',[1 z]);
if nii_out, print('dhb.eps','-depsc2','-r300'), end

%% Output parameter niftis
if nii_out
    
    save_avw(r2p, 'r2p', 'f', scales)
    save_avw(dbv, 'dbv', 'f', scales)
    save_avw(oef, 'oef', 'f', scales)
    save_avw(dhb, 'dhb', 'f', scales)

    saveas(imgH_r2p, 'r2p.fig')
    saveas(imgH_dbv, 'dbv.fig')
    saveas(imgH_oef, 'oef.fig')
    saveas(imgH_dhb, 'dhb.fig')
   
end

end
end
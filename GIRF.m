%% Do Generalized IRF
clear
load main.mat
n_shocks = 2; % unc and tfp shock
irfperiod = 40;

%% Get draws from previous simulations
n_worlds = 1000;
startpoints = cell(n_worlds,1);
innov = cell(n_worlds,1);
rng('default');
rng(2015);
for i_world = 1:n_worlds
	startpoints{i_world}.K = Ksim(end-i_world+1);
	startpoints{i_world}.ssigmaxind = states_sim.ssigmaxind(end-i_world+1);
	startpoints{i_world}.zind = states_sim.zind(end-i_world+1);
	startpoints{i_world}.dist_k = dist_k(:,:,end-i_world+1);
	innov{i_world}.rand_z = rand(1,T);
	innov{i_world}.rand_unc = rand(1,T);
end

%% Packing parameters
package.K_grid = K_grid;
package.fine_grid = fine_grid;
package.noinvest_ind_fine = noinvest_ind_fine;
package.ssigmax_grid = ssigmax ; % careful here name difference
package.markup_grid = markup_grid;
package.z_grid = Z ;
package.q_grid = q_grid ;
package.pphi_c = pphi_C ;
package.pphi_K = pphi_K ;
package.pphi_q = pphi_q ;
package.pphi_ttheta = pphi_ttheta ;
package.pphi_tthetaq = pphi_tthetaq ;
package.PX_low = PX_low ;
package.PX_high = PX_high ;
package.X = X;
package.nz = nz ;
package.nx = nx ;
package.ssigmax_cdf = ssigmax_cdf;
package.z_cdf = z_cdf;

%% Compute GIRFs
controls_unc_high   = cell(n_worlds,1);
controls_unc_low    = cell(n_worlds,1);
states_unc_high     = cell(n_worlds,1);
states_unc_low      = cell(n_worlds,1);
parfor i_world = 1:n_worlds
	[states_unc_high{i_world},controls_unc_high{i_world}] = simulation(irfperiod,startpoints{i_world},innov{i_world},'unc_high',kopt_fine,active_fine,package);
	[states_unc_low{i_world} ,controls_unc_low{i_world} ] = simulation(irfperiod,startpoints{i_world},innov{i_world},'unc_low',kopt_fine,active_fine,package);
	i_world
end
save GIRF.mat

%% Plotting for uncertainty shock
q_panel_unc_high      = zeros(n_worlds,irfperiod);
q_panel_unc_low       = zeros(n_worlds,irfperiod);
C_panel_unc_high      = zeros(n_worlds,irfperiod);
C_panel_unc_low       = zeros(n_worlds,irfperiod);
w_panel_unc_high      = zeros(n_worlds,irfperiod);
w_panel_unc_low       = zeros(n_worlds,irfperiod);
ttheta_panel_unc_high = zeros(n_worlds,irfperiod);
ttheta_panel_unc_low  = zeros(n_worlds,irfperiod);
inv_panel_unc_high    = zeros(n_worlds,irfperiod);
inv_panel_unc_low     = zeros(n_worlds,irfperiod);
GDP_panel_unc_high    = zeros(n_worlds,irfperiod);
GDP_panel_unc_low     = zeros(n_worlds,irfperiod);
parfor i_world = 1:n_worlds
	q_panel_unc_high(i_world,:) = controls_unc_high{i_world}.q;
	q_panel_unc_low(i_world,:)  = controls_unc_low{i_world}.q ;
	C_panel_unc_high(i_world,:) = controls_unc_high{i_world}.C;
	C_panel_unc_low(i_world,:)  = controls_unc_low{i_world}.C;
	w_panel_unc_high(i_world,:) = controls_unc_high{i_world}.w;
	w_panel_unc_low(i_world,:)  = controls_unc_low{i_world}.w;
	ttheta_panel_unc_high(i_world,:) = controls_unc_high{i_world}.ttheta;
	ttheta_panel_unc_low(i_world,:)  = controls_unc_low{i_world}.ttheta;
	inv_panel_unc_high(i_world,:) = controls_unc_high{i_world}.inv;
	inv_panel_unc_low(i_world,:)  = controls_unc_low{i_world}.inv;
	GDP_panel_unc_high(i_world,:) = controls_unc_high{i_world}.GDP;
	GDP_panel_unc_low(i_world,:)  = controls_unc_low{i_world}.GDP;
end
h1 = figure;
subplot(3,2,1);
plot(1:irfperiod,mean(q_panel_unc_high)-mean(q_panel_unc_low));
title('Price')
subplot(3,2,2);
plot(1:irfperiod,mean(C_panel_unc_high)-mean(C_panel_unc_low));
title('Consumption')
subplot(3,2,3);
plot(1:irfperiod,mean(w_panel_unc_high)-mean(w_panel_unc_low));
title('Wage')
subplot(3,2,4);
plot(1:irfperiod,mean(ttheta_panel_unc_high)-mean(ttheta_panel_unc_low));
title('Tightness')
subplot(3,2,5);
plot(1:irfperiod,mean(inv_panel_unc_high)-mean(inv_panel_unc_low));
title('Investment')
subplot(3,2,6);
plot(1:irfperiod,mean(GDP_panel_unc_high)-mean(GDP_panel_unc_low));
title('GDP')
print(h1,'GIRF_unc.eps','-depsc2')

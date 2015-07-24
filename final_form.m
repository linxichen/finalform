%% Changelog
% 05/05/2015 Linxi: Add intensive margin to simplified version. This should
% be the only version we use from now. Heavily comment this one.

% 05/19/2015 General equilibrium. Use quarterly calibration.
% Not much investment frictions yet.
% Key assumptions:
% 1. resale price is 0. Means irreversible capital. DONT CHANGE THIS!
% otherwise have to rewrite a lot of codes
% 2. wage = psi_n*c the wage supply equation
% 3. General eqm, CARA
% 4. Uses steady state calibration file

% 6/29 quadratic adjustment cost
%      k = (1-delta)k+I-eta*(k*(I/K)^2) (Bloom 2009, p13)

% 07/21 eliminated interpolation
% 07/23 major overhual of simulation part

%% Housekeeping
clear; close all; clc;
deep_para_quarterly;
diary('log.txt');

%% Accuracy control
nk = 100; % number of grid points on capital stock
nfine = 100; % number of grid points for policy functions and simulation
nx = 7; % number of grid points on idiosyncractic prod.
nz = 7; % number of grid points on aggregate productivity
ns = nx*nz*2; % total dimension of exo state, (idio, agg, ssigmax)
nK = 15; % agg capital state
nq = 15; % number of price points inv producer can set
nmarkup = 15;
m = 2.5; % support is m s.d. away from mean
tol = 1e-3; % when to stop VFI and KS
outer_tol = 1e-3;
maxiter = 100;
damp = 0.5; % help with outer convergence
T = 3000; % How long to simulate
burnin = ceil(0.1*T); % how many periods to discard to get rid of dependence on initial state

%% Grids and stuff derived immediately from parameters
[Z,PZ] = tauchen(nz,0,rrhoz,ssigmaz,m); Z = exp(Z); % agg TFP
ssigmax = [ssigmax_low,ssigmax_high];
[X,PX_high] = tauchen(nx,0,rrhox,ssigmax_high,m);
%=========================================================================%
% This PX_high(i,j) gives the approximate prob of P(X_t+1(j) | X_t = X(i))
% Follows X_t+1 = rrhox*X(i) + ssigma_high* N(0,1)
%=========================================================================%
PX_low = tauchen_givengrid(0,rrhox,ssigmax_low,X); X = exp(X);
P = zeros(nx*nz*2,nx*nz*2);  % The transition of exo state vec, all independent
low = 1;
high = 2;
% Below find the transition prob of three indie variable
for i = 1:ns
    [i_x,i_z,i_ssigmax] = ind2sub([nx nz 2],i);
    for j = 1:ns
        [j_x,j_z,j_ssigmax] = ind2sub([nx nz 2],j);
        P(i,j) = PZ(i_z,j_z)*Pssigmax(i_ssigmax,j_ssigmax)*( (j_ssigmax==low)*PX_low(i_x,j_x)+ (j_ssigmax==high)*PX_high(i_x,j_x) );
    end
end
ssigmax_cdf = cumsum(Pssigmax,2);
z_cdf = cumsum(PZ,2);

%=========================================================================%
% inv price grid. Since MC is 1, the price must be larger than 1
% From the two period model, q = (xxi - bbeta*(1-ddelta))/(xxi - 1),
% The tail index has to be larger than 1/(1-aalpha) = 1.4286, this implies
% q = 1.05.
%=========================================================================%
markupmin = 1.00;
markupmax = 1.2;
markup_grid = linspace(markupmin,markupmax,nmarkup);
q_grid = linspace(0.4,1.8,nq);
% Capital grid
min_k = 1;
max_k = 10;
k_grid = zeros(nk,1);
for i_k = 1:nk
    k_grid(i_k) = max_k*(1-ddelta)^(nk-i_k);
end
% k_grid = linspace(min_k,max_k,nk)';
% fine_grid = linspace(k_grid(1),k_grid(nk),nfine)';
fine_grid = k_grid;
noinvest_ind = ones(nk,1); % for each k, the index of tmr k if no invest
for i_k = 1:nk
    [~,noinvest_ind(i_k)] = min(abs(k_grid-(1-ddelta)*k_grid(i_k)));
end
noinvest_ind_fine = ones(nfine,1); % for each k, the index of tmr k if no invest
for i_k = 1:nfine
    [~,noinvest_ind_fine(i_k)] = min(abs(fine_grid-(1-ddelta)*fine_grid(i_k)));
end
inv_mat = repmat(k_grid',nk,1)-(1-ddelta)*repmat(k_grid,1,nk);
%=========================================================================%
% This inv_mat(i,j) gives inv amount needed to reach k_t+1 = k(j) when
% k_t = k(i)
%=========================================================================%
pos_inv = inv_mat>0;
neg_inv = inv_mat<=0;
inv_mat_fine = repmat(fine_grid',nfine,1)-(1-ddelta)*repmat(fine_grid,1,nfine);
pos_inv_fine = inv_mat_fine>0;
neg_inv_fine = inv_mat_fine<=0;
K_grid = linspace(k_grid(1),k_grid(nk),nK)'; % Aggregate capital grid

%% Aggregate rules and its parameters
% Rule of q as function of aggregate states
pphi_qC = log(mean(q_grid)); % constant term
pphi_qK = 0; % w.r.t agg K
pphi_qz = 0; % w.r.t agg TFP
pphi_qssigmax = 0; % w.r.t uncertainty
pphi_q = [pphi_qC,pphi_qK,pphi_qssigmax,pphi_qz];

pphi_KK = 9.870171e-01;
pphi_KC = log(mean(k_grid));
pphi_Kz = 0;
pphi_Kssigmax = 0;% Aggregate Law of motion for aggregate capital
pphi_K = [pphi_KC,pphi_KK,pphi_Kssigmax,pphi_Kz];

pphi_CC = log(1.3);
pphi_CK = 0.0;
pphi_Cz = 0.0;
pphi_Cssigmax = 0.0;
pphi_C = [pphi_CC,pphi_CK,pphi_Cz,pphi_Cssigmax];

pphi_tthetaC = log(1); % tightness ratio depends on q
pphi_tthetaK = 0.0;
pphi_tthetaz = 0.0;
pphi_tthetassigmax = 0.0;
pphi_ttheta = [pphi_tthetaC,pphi_tthetaK,pphi_tthetaz,pphi_tthetassigmax];
pphi_tthetaq = 0.1;% lower q -- more firms invest -- lower ttheta
if (exist('aggrules.mat','file') == 2)
    load aggrules.mat;
end

%% Initialize value functions
W_old = ones(nk,ns,nK,nq); % value of matching with investment goods producer after paying the search cost
W_new = W_old;
U_old = ones(nk,ns,nK); % value of not going to search, not related to current q
U_new = U_old;
V_old = ones(nk,ns,nK,nq); % maximized value after discrete choice
V_new = ones(nk,ns,nK,nq); %
if (exist('valuefunctions.mat','file') == 2)
    load valuefunctions.mat
else
    W_old = ones(nk,ns,nK,nq); % value of matching with investment goods producer after paying the search cost
    W_new = W_old;
    U_old = ones(nk,ns,nK); % value of not going to search, not related to current q
    U_new = U_old;
    V_old = ones(nk,ns,nK,nq); % maximized value after discrete choice
    V_new = ones(nk,ns,nK,nq); %
end
if isequal(size(W_old),[nk,ns,nK,nq])
else
    W_old = ones(nk,ns,nK,nq); % value of matching with investment goods producer after paying the search cost
    W_new = W_old;
    U_old = ones(nk,ns,nK); % value of not going to search, not related to current q
    U_new = U_old;
    V_old = ones(nk,ns,nK,nq); % maximized value after discrete choice
    V_new = ones(nk,ns,nK,nq); %
end

% profit_fine = zeros(nfine,ns,nq);
V_new_fine = zeros(nfine,ns,nK,nq);
EV_new_fine = zeros(nfine,ns,nK,nq);
W_new_fine = zeros(nfine,ns,nK,nq);
U_new_fine = zeros(nfine,ns,nK);
koptind_active = zeros(nk,ns,nK,nq);
koptind = zeros(nk,ns,nK,nq);
kopt_active_fine = zeros(nfine,ns,nK,nq);
active = zeros(nk,ns,nK,nq);

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

% New Discrete Simulation? simulate the innovation terms
x_cdf_low = cumsum(PX_low,2);
x_cdf_high = cumsum(PX_high,2);
rng('default');
rng(2015);
innov.rand_z = rand(1,T);
innov.rand_unc = rand(1,T);
startpoints.K = mean(K_grid);
startpoints.ssigmaxind = low;
startpoints.zind = ceil(nz/2);
startpoints.dist_k = zeros(nfine,nx);
startpoints.dist_k(ceil(nfine/2),:) = ones(1,nx)/nx;

tic
%% Main Body of KS iter
outer_diff = 10;
outer_iter = 0;

while ((outer_diff > outer_tol) && (outer_iter < maxiter))
	%============ VFI Begins==============================================%
	% Inner loop
	err = 10;
	iter = 0;
	while err > tol
		for i_s = 1:ns
			[i_z,i_x,i_ssigmax] = ind2sub([nz nx 2],i_s);
			for i_K = 1:nK
				% What is agg state today
				log_aggstate = [1; log(K_grid(i_K)); log(ssigmax(i_ssigmax)); log(Z(i_z))];
				% Forecast future values
				C = exp(pphi_C*log_aggstate);
				w = ppsi_n*C;
				[~,i_Kplus] = min(abs(K_grid-exp(pphi_K*log_aggstate)));
				[~,i_qplus] = min(abs(q_grid-exp(pphi_q*log_aggstate)));
				%=========Take as given wage,find profit in c units=======%
				L = (w*k_grid.^(-aalpha)/Z(i_z)/X(i_x)/v).^(1/(v-1));
				profit = Z(i_z)*X(i_x)*k_grid.^aalpha.*L.^v - w.*L;
				L = repmat(L,1,nk);
				%=========================================================$

				% Find expected value across i_splus
				EV_noinvest = V_old(noinvest_ind,:,i_Kplus,i_qplus)*P(i_s,:)'; % EV_oninvest(i) gives expected V given k_t+1 = K(i), optimally!

				% U has no choice, so no need to find max in a loop.
				llambda = C^(-1);
				U_new(:,i_s,i_K) = llambda*profit + bbeta*EV_noinvest;

				% Find W
				EV_invest = V_old(:,:,i_Kplus,i_qplus)*P(i_s,:)';
				%=========================================================%
				% This EV_invest(i) is expected value when kplus=K(i)
				% Not optimally! We need to find maximizing kplus
				%=========================================================%
				for i_q = 1:nq
					eeta = 0.5;
					convexadj = eeta*repmat(k_grid,1,nk).*(inv_mat./repmat(k_grid,1,nk)).^2;
					ttheta = exp(pphi_ttheta*log_aggstate+pphi_tthetaq*log(q_grid(i_q))); % tightness for this q
					% mmu = 1./(1+tttheta_old.^(-aalpha0)).^(1/aalpha0);
					mmu = aalpha0*ttheta^aalpha1;
					% mmu = 1./((1+ttheta.^(aalpha0)).^(1/aalpha0));
					netprofit = repmat(profit,1,nk)-mmu*q_grid(i_q)*(inv_mat).*(pos_inv+neg_inv*pphi)-mmu*convexadj;
					candidate = llambda*netprofit  + bbeta*(mmu*repmat(EV_invest',nk,1)+(1-mmu)*repmat(EV_noinvest,1,nk));
					%=====================================================%
					% candidate(i_k,i_kplus) is the thing inside max
					% operator. expected value tomorrow depends only on
					% capital tomorrow and aggregate state, therefore I use
					% repmat here. netprofit doesn't depend on value
					% function so it is pre-computed outside for loops
					%=====================================================%
					[W_new(:,i_s,i_K,i_q),koptind_active(:,i_s,i_K,i_q)] = max(candidate,[],2);
					% W_new(:,i_s,i_K,i_q) = mmu*W_temp+(1-mmu)*U_new(:,i_s,i_K); % match with prob. mmu
					V_new(:,i_s,i_K,i_q) = max(W_new(:,i_s,i_K,i_q),U_new(:,i_s,i_K));
				end
			end
		end

		err = norm([V_old(:);W_old(:);U_old(:)]-[V_new(:);W_new(:);U_new(:)],Inf);
		V_old = V_new;
		W_old = W_new;
		U_old = U_new;
		iter = iter + 1;
		if mod(iter,1) == 0
			disp_text = sprintf('KS Iter = %d, KS err = %d, Current VFI Iter = %d, err = %d',outer_iter,outer_diff,iter,err);
			disp(disp_text);
		end
	end
	%============ VFI Ends================================================%

	active = W_new > repmat(U_new,1,1,1,nq);
	koptind = repmat(noinvest_ind,1,ns,nK,nq).*(1-active) + active.*koptind_active;
	kopt_active = k_grid(koptind_active);
	kopt = k_grid(koptind);
	plot(k_grid,kopt(:,sub2ind([nx nz 2],ceil(nx/2),ceil(nz/2),2),ceil(nK/2),ceil(nq/2))-(1-ddelta)*k_grid)
	save('valuefunctions.mat','V_new','W_new','U_new','V_old','W_old','U_old')

	%     % Interpolate on finer grid
	%     for i_q = 1:nq
	%         for i_K = 1:nK
	%             for i_s = 1:ns
	%                 U_new_fine(:,i_s,i_K) = interp1(k_grid,U_new(:,i_s,i_K),fine_grid,'linear')';
	%                 W_new_fine(:,i_s,i_K,i_q) = interp1(k_grid,W_new(:,i_s,i_K,i_q),fine_grid,'linear')';
	%                 kopt_active_fine(:,i_s,i_K,i_q) = interp1(k_grid,kopt_active(:,i_s,i_K,i_q),fine_grid,'linear')';
	%             end
	%         end
	%     end
	U_new_fine = U_new;
	W_new_fine = W_new;
	kopt_active_fine = kopt_active;
	active_fine = W_new_fine > repmat(U_new_fine,1,1,1,nq);
	kopt_fine = (1-ddelta)*repmat(fine_grid,1,ns,nK,nq).*(1-active_fine) + active_fine.*kopt_active_fine;
	plot(fine_grid,kopt_fine(:,sub2ind([nx nz 2],ceil(nx/2),ceil(nz/2),2),ceil(nK/2),1)-(1-ddelta)*fine_grid)

	%% Given individual policies, simulate a large panel to update aggregate law of motion
	[states_sim,controls_sim] = simulation(T,startpoints,innov,'allshock',kopt_fine,active_fine,package);
	Ksim = states_sim.K;
	ssigmaxsim = states_sim.ssigmax;
	zsim = states_sim.z;
	qsim = controls_sim.q;
	tthetasim = controls_sim.ttheta;
	Csim = controls_sim.C;

	%% Regress to get coefficients of K law
	XX = [ones(T-burnin-1,1) log(Ksim(burnin+1:T-1))' log(ssigmaxsim(burnin+1:T-1))' log(zsim(burnin+1:T-1))'];
	Y = log(Ksim(2+burnin:T)');
	bbeta_K = (XX'*XX)\(XX'*Y);
	e = Y-XX*bbeta_K;
	ytilde = Y-mean(Y);
	Rsq_K = 1-(e'*e)/(ytilde'*ytilde);
	pphi_K_new = damp*bbeta_K' + (1-damp)*pphi_K;

	% Regress to get q law
	Y = log(qsim(1+burnin:T-1))';
	bbeta_q = (XX'*XX)\(XX'*Y);
	e = Y-XX*bbeta_q;
	ytilde = Y-mean(Y);
	Rsq_q = 1-(e'*e)/(ytilde'*ytilde);
	pphi_q_new = damp*bbeta_q' + (1-damp)*pphi_q;

	% Regress to get C law
	Y = log(Csim(1+burnin:T-1))';
	bbeta_C = (XX'*XX)\(XX'*Y);
	e = Y-XX*bbeta_C;
	ytilde = Y-mean(Y);
	Rsq_C = 1-(e'*e)/(ytilde'*ytilde);
	pphi_C_new = damp*bbeta_C' + (1-damp)*pphi_C;

	% Regress to get ttheta law
	Y = log(tthetasim(1+burnin:T-1))';
	XX1 = [XX,log(qsim(1+burnin:T-1))'];
	bbeta_ttheta = (XX1'*XX1)\(XX1'*Y);
	e = Y-XX1*bbeta_ttheta;
	ytilde = Y-mean(Y);
	Rsq_ttheta = 1-(e'*e)/(ytilde'*ytilde);
	pphi_ttheta_new = damp*bbeta_ttheta(1:end-1)' + (1-damp)*pphi_ttheta;
	pphi_tthetaq_new = damp*bbeta_ttheta(end) + (1-damp)*pphi_tthetaq;

	% Update mmu_old as well
	outer_diff = norm([pphi_K,pphi_q,pphi_C,pphi_ttheta,pphi_tthetaq]-[pphi_K_new,pphi_q_new,pphi_C_new,pphi_ttheta_new,pphi_tthetaq_new],Inf);

	% Update mmu_old as well
	pphi_K = pphi_K_new;
	pphi_q = pphi_q_new;
	pphi_C = pphi_C_new;
	pphi_ttheta = pphi_ttheta_new;
	pphi_tthetaq = pphi_tthetaq_new;

	%% Plot something
	%     i_x = 4;
	%     i_K = 3;
	%     i_q = 3;
	%     figure
	%     plot(k_grid,W_old(:,i_x,i_K,i_q),'-r',k_grid,U_old(:,i_x,i_K),'-b');
	%     figure
	%     plot(k_grid,V_old(:,i_x,i_K,i_q));
	%     figure
	%     plot(k_grid,kopt(:,i_x,i_K,i_q));
	%     figure
	%     plot(q_grid,revenue(:,T),'-.r',q_grid,revenue(:,T-7),'b');
	%     figure
	%     plot(q_grid,demand(:,T),'-.r',q_grid,demand(:,T-7),'b');

	%%
	outer_iter = outer_iter + 1;
	disp_text = sprintf('Rsq_K = %d, Rsq_q = %d,Rsq_C = %d,Rsq_ttheta = %d',Rsq_K,Rsq_q,Rsq_C,Rsq_ttheta);
	disp(disp_text);
	disp_text = sprintf('log(q) = %d + %d * log(K) + %d * log(ssigmax)+%d * log(z)',pphi_q);
	disp(disp_text);
	disp_text = sprintf('log(Kplus) = %d + %d * log(K) + %d * log(ssigmax)+%d * log(z)',pphi_K);
	disp(disp_text);
	disp_text = sprintf('log(C) = %d + %d * log(K) + %d * log(ssigmax)+%d * log(z)',pphi_C);
	disp(disp_text);
	disp_text = sprintf('log(ttheta) = %d + %d * log(K) + %d * log(ssigmax)+%d * log(z)+%d*log(q)',pphi_C,pphi_tthetaq);
	disp(disp_text);
	disp_text = sprintf('KS Iter = %d, KS err = %d, Current VFI Iter = %d, err = %d',outer_iter,outer_diff,iter,err);
	disp(disp_text);
	disp('===============================');
	save('aggrules.mat','pphi_K','pphi_q','pphi_C','pphi_ttheta');
	save('Rsq.mat','Rsq_K','Rsq_q','Rsq_C','Rsq_ttheta');

end


toc



%% Decompose Demand Curve
% According to policy functions, find the optimal q
% Fix the distribution of firms the last period T
% t = T;
% [~,i_K] = min(abs((Ksim(t)-K_grid)));
% revenue_lowtfp = zeros(1,nq);
% demand_lowtfp = revenue_lowtfp;
% whichs = zeros(1,nx);
% for i_x = 1:nx
%     whichs(i_x) = sub2ind([nz nx 2],1,i_x,ssigmaxsim(t));
% end
% for i_q = 1:nq
%     tot_profit_grid(:,:,i_q) = (q_grid(i_q)-psi)*dist_k(:,:,t).*((0.02+ddelta)*repmat(fine_grid,1,nx)).*(active_fine(:,whichs,i_K,i_q));
%     % tot_revenue_grid(tot_revenue_grid<0) = 0;
%     demand_grid = dist_k(:,:,t).*((0.02+ddelta)*repmat(fine_grid,1,nx)).*(active_fine(:,whichs,i_K,i_q));
%     revenue_lowtfp(i_q) = sum(vec(tot_profit_grid(:,:,i_q)));
%     demand_lowtfp(i_q) = sum(demand_grid(:));
% end
% figure
% mesh(X,fine_grid,tot_profit_grid(:,:,1))
%
% t = T;
% [~,i_K] = min(abs((Ksim(t)-K_grid)));
% revenue_hightfp = zeros(1,nq);
% demand_hightfp = revenue_hightfp;
% whichs = zeros(1,nx);
% for i_x = 1:nx
%     whichs(i_x) = sub2ind([nz nx 2],nz,i_x,ssigmaxsim(t));
% end
% for i_q = 1:nq
%     tot_profit_grid(:,:,i_q) = (q_grid(i_q)-psi)*dist_k(:,:,t).*((0.02+ddelta)*repmat(fine_grid,1,nx)).*(active_fine(:,whichs,i_K,i_q));
%     % tot_revenue_grid(tot_revenue_grid<0) = 0;
%     demand_grid = dist_k(:,:,t).*((0.02+ddelta)*repmat(fine_grid,1,nx)).*(active_fine(:,whichs,i_K,i_q));
%     revenue_hightfp(i_q) = sum(vec(tot_profit_grid(:,:,i_q)));
%     demand_hightfp(i_q) = sum(demand_grid(:));
% end
%

figure
mesh(X,fine_grid,tot_profit_grid(:,:,1))

% plot(q_grid,revenue_lowtfp,'b',q_grid,revenue_hightfp,'r')
save main.mat

% checkresults;

%% Understand results
period = T-2;

figure_dist = figure;
surf(X,fine_grid,real(dist_k(:,:,end-1)))
savefig(figure_dist,'distribution.fig');

figure_activemeasure = figure;
whichs = zeros(1,nx);
for i_x = 1:nx
	whichs(i_x) = sub2ind([nz nx 2],zindsim(period),i_x,ssigmaxsim(period));
end
surf(X,fine_grid,real(dist_k(:,:,end).*(active_fine(:,whichs,i_K,i_q))))
savefig(figure_activemeasure,'active_measure.fig');

figure_ksim = figure;
plot(1:T,Ksim,1:T,qsim);
xlabel('Time');
ylabel('Aggregate Capital');
savefig(figure_ksim,'Ksim.fig');
legend('Capital','Inv Good Price');

figure_csim = figure;
plot(1:T,Csim,1:T,ssigmaxsim)
xlabel('Time');
legend('Consumption','Uncertainty');
savefig(figure_csim,'csim.fig');

figure_tthetasim = figure;
plot(1:T,qsim,1:T,tthetasim)
xlabel('Time');
legend('Inv Price','Tightness');
savefig(figure_tthetasim,'tthetasim.fig');

figure_invsim = figure;
plot(2:T,Ksim(2:T)-(1-ddelta)*Ksim(1:T-1),2:T,qsim(2:T))
xlabel('Time');
legend('Inv Quantity','Tightness');
savefig(figure_tthetasim,'tthetasim.fig');

mmusim = aalpha0*real(tthetasim).^(aalpha1);
activesim = 1./tthetasim;

invsim = real(Ksim(2:T)-(1-ddelta)*Ksim(1:T-1));
inv_mean = real(mean(invsim));
[~,inv_cyc] = hpfilter(invsim./inv_mean-1,1600);
std(inv_cyc)
diary off;

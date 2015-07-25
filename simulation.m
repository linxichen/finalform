function [states_sim,controls_sim] = simulation(T,startpoints,innov,whichshock,kopt_fine,active_fine,package)
% Simulate the economy given aggregate rules and invidual policies.
% state variables are in the order of [K,dist_k,z,ssigmax]
% control variables are in the order of [C,w,q,]
% Input:
%   T: periods of simulation.
%   startpoints: a structure that contains the starting state vector. e.g. startpoints.K is the initial
%                aggregate capital level
%   innov: nshock-by-T matrix of innovation terms. Usually N(0,1) or U(0,1)
%   whickshock: string specifies which shock(s) will hit the economy at t=1. All other periods are random.
% Output:
%   states_sim: simulated sequence of state variables. Stored in an array of structure.
%   controls_sim: simulated sequence of control variables. Stored in an array of structure.

%% Read parameters (things that don't change at runtime) and unpack stuff
deep_para_quarterly;
K_grid                   = package.K_grid;
fine_grid                = package.fine_grid;
noinvest_ind_fine        = package.noinvest_ind_fine;
ssigmax_grid             = package.ssigmax_grid; % careful here name difference
markup_grid              = package.markup_grid; % careful here name difference
z_grid                   = package.z_grid;
q_grid                   = package.q_grid;
pphi_C                   = package.pphi_c;
pphi_ttheta              = package.pphi_ttheta;
pphi_tthetaq             = package.pphi_tthetaq;
PX_low                   = package.PX_low;
PX_high                  = package.PX_high;
X                        = package.X;
nz                       = package.nz;
nx                       = package.nx;
z_cdf                    = package.z_cdf;
ssigmax_cdf              = package.ssigmax_cdf;

rand_z                   = innov.rand_z;
rand_unc                 = innov.rand_unc;

%% Who's who
low = 1;
high = 2;

%% Preallocation state and control simulation arrays
nfine = length(fine_grid);
nmarkup = length(markup_grid);
K_sim                    = zeros(1,T);
dist_k_sim               = zeros(nfine,nx,T);
zind_sim                 = zeros(1,T);
ssigmaxind_sim           = zeros(1,T);

C_sim                    = zeros(1,T);
w_sim                    = zeros(1,T);
q_sim                    = zeros(1,T);
ttheta_sim               = zeros(1,T);
revenue_sim              = zeros(nmarkup,T);
demand_sim               = zeros(nmarkup,T);

%% Assign Starting states
K_sim(1) = startpoints.K;
dist_k_sim(:,:,1) = startpoints.dist_k;
zind_sim(1) = startpoints.zind;
ssigmaxind_sim(1) = startpoints.ssigmaxind;

%% Simulate forward
for t = 1:T
	% Find aggregate stuff today, assuming all state variables are assigned
	K_sim(t) = sum(vec(dist_k_sim(:,:,t).*repmat(fine_grid,1,nx)));
	[~,i_K] = min(abs(K_sim(t)-K_grid));
	% Given uncertainty and agg TFP, find the agg exo state index of each indiviual
	% prod. level. He uses this to forecast expected value function
	whichs = zeros(1,nx);
	for i_x = 1:nx
		whichs(i_x) = sub2ind([nx nz 2],i_x,zind_sim(t),ssigmaxind_sim(t)); % previously in the wrong order!
	end
	z       = z_grid(zind_sim(t));
	ssigmax = ssigmax_grid(ssigmaxind_sim(t));

	% Find transition matrix according to today's state
	if ssigmaxind_sim(t) == low
		whichprob = PX_low;
	elseif ssigmaxind_sim(t) == high
		whichprob = PX_high;
	end

	% Find consumption and wage today
	log_aggstate = [1; log(K_sim(t)); log(ssigmax); log(z)];
	C = exp(pphi_C*log_aggstate);
	w = ppsi_n*C;

	% According to policy functions, find the optimal q
	for i_markup = 1:nmarkup
		[~,i_q] = min(abs(q_grid-markup_grid(i_markup)*w));
		ttheta_temp = exp(pphi_ttheta*log_aggstate+pphi_tthetaq*log(q_grid(i_q)));% tightness ratio given q and states
		mmu_temp = aalpha0*ttheta_temp^aalpha1;
		if (mmu_temp>1)
			warning('mmu > 1 encountered.')
		end
		tot_profit_grid = mmu_temp*(markup_grid(i_markup)*w-w)*dist_k_sim(:,:,t).*(kopt_fine(:,whichs,i_K,i_q)-(1-ddelta)*repmat(fine_grid,1,nx)).*(active_fine(:,whichs,i_K,i_q));
		demand_grid = mmu_temp*dist_k_sim(:,:,t).*(kopt_fine(:,whichs,i_K,i_q)-(1-ddelta)*repmat(fine_grid,1,nx)).*(active_fine(:,whichs,i_K,i_q));
		revenue_sim(i_markup,t) = sum(tot_profit_grid(:));
		demand_sim(i_markup,t) = sum(demand_grid(:));
	end
	[~,i_markupmax] = max(revenue_sim(:,t));
	qmax = w*markup_grid(i_markupmax);
	[~,i_qmax] = min(abs(q_grid-qmax));
	q_sim(t) = qmax;

	% Acquire q from rule
	% qmax = exp(pphi_q*log_aggstate);
	% [~,i_qmax] = min(abs(q_grid-qmax));

	% Evolution under the argmax q
	if (t<=T-1)
		output = 0;
		ttheta_temp = exp(pphi_ttheta*log_aggstate+pphi_tthetaq*log(qmax));% tightness ratio given q and states
		mmu_temp = aalpha0*ttheta_temp^aalpha1;
		if mmu_temp > 1
			warning('mmu > 1 encountered. Debug.');
		end
		i_z = zind_sim(t);
		for i_k = 1:nfine
			for i_x = 1:nx
				%======GE: Find output on each state==================%
				L = (w*fine_grid(i_k).^(-aalpha)/z/X(i_x)/v).^(1/(v-1));
				output = output + dist_k_sim(i_k,i_x,t)*(z*X(i_x)*fine_grid(i_k).^aalpha.*L.^v); % previously i_z is not correctly created from zsim(t). It stucked at VFI stage!!!
				%======GE: Find output on each state==================%

				i_s = sub2ind([nx nz 2],i_x,i_z,ssigmaxind_sim(t));
				kplus = kopt_fine(i_k,i_s,i_K,i_qmax);

				% Assign mass to tomorrow's distribution
				if active_fine(i_k,i_x,i_K,i_qmax) == 1 % previously forgot the i_K,i_qmax part
					if (kplus>=fine_grid(1) && kplus<fine_grid(end))
						lower_ind = find(fine_grid<=kplus,1,'last');
						upper_ind = lower_ind + 1;
						denom = fine_grid(upper_ind)-fine_grid(lower_ind);
						for i_xplus = 1:nx
							dist_k_sim(lower_ind,i_xplus,t+1) = dist_k_sim(lower_ind,i_xplus,t+1) + mmu_temp*whichprob(i_x,i_xplus)*dist_k_sim(i_k,i_x,t)*(fine_grid(upper_ind)-kplus)/denom;
							dist_k_sim(upper_ind,i_xplus,t+1) = dist_k_sim(upper_ind,i_xplus,t+1) + mmu_temp*whichprob(i_x,i_xplus)*dist_k_sim(i_k,i_x,t)*(kplus-fine_grid(lower_ind))/denom;
							dist_k_sim(noinvest_ind_fine(i_k),i_xplus,t+1) = dist_k_sim(noinvest_ind_fine(i_k),i_xplus,t+1)+(1-mmu_temp)*whichprob(i_x,i_xplus)*dist_k_sim(i_k,i_x,t);
						end
					elseif (kplus<fine_grid(1))
						for i_xplus = 1:nx
							dist_k_sim(1,i_xplus,t+1) = dist_k_sim(1,i_xplus,t+1) + mmu_temp*whichprob(i_x,i_xplus)*dist_k_sim(i_k,i_x,t);
							dist_k_sim(noinvest_ind_fine(i_k),i_xplus,t+1) = dist_k_sim(noinvest_ind_fine(i_k),i_xplus,t+1)+(1-mmu_temp)*whichprob(i_x,i_xplus)*dist_k_sim(i_k,i_x,t);
						end
					elseif (kplus>=fine_grid(nfine))
						for i_xplus = 1:nx
							dist_k_sim(nfine,i_xplus,t+1) = dist_k_sim(nfine,i_xplus,t+1) + mmu_temp*whichprob(i_x,i_xplus)*dist_k_sim(i_k,i_x,t);
							dist_k_sim(noinvest_ind_fine(i_k),i_xplus,t+1) = dist_k_sim(noinvest_ind_fine(i_k),i_xplus,t+1)+(1-mmu_temp)*whichprob(i_x,i_xplus)*dist_k_sim(i_k,i_x,t);
						end
					end
				else
					for i_xplus = 1:nx
						dist_k_sim(noinvest_ind_fine(i_k),i_xplus,t+1) = dist_k_sim(noinvest_ind_fine(i_k),i_xplus,t+1)+whichprob(i_x,i_xplus)*dist_k_sim(i_k,i_x,t);
					end
				end
			end
		end

		% Check if distribution makes sense
		if min(vec(dist_k_sim(:,:,t+1))) < 0
			warning('somewhere theres negative mass')
		end
		if sum(vec(dist_k_sim(:,:,t+1))) ~= 1
			sprintf('Mass sums up to %d',sum(vec(dist_k_sim(:,:,t+1))));
		end
	end

	% Eventually find implied consumption
	C_sim(t) = output;
	% also find implied tightness ratio
	measure_active = sum(sum(dist_k_sim(:,:,t).*active_fine(:,whichs,i_K,i_qmax)));% how many firms choose to search
	if measure_active == 0
		ttheta_sim(t) = 999999;
	else
		ttheta_sim(t) = 1/measure_active;
	end
	q_sim(t) = qmax;
	w_sim(t) = C_sim(t)*ppsi_n;

	if t == 1 % apply irf shock only in the first period (t==2)
		if strcmp(whichshock,'tfp')
			z_temp                      = exp(rrhoz*log(z_grid(zind_sim(t)))+ssigmaz*1);
			[~,i_ztemp]                 = min(abs(z_grid-z_temp));
			zind_sim(t+1)               = i_ztemp;
			ssigmaxind_sim(t+1)         = markov_draw(ssigmaxind_sim(t),ssigmax_cdf,rand_unc(t));
		elseif strcmp(whichshock,'unc_high')
			ssigmaxind_sim(t+1)         = high;
			zind_sim(t+1)               = markov_draw(zind_sim(t)      ,z_cdf      ,rand_z(t));
		elseif strcmp(whichshock,'unc_low')
			ssigmaxind_sim(t+1)         = low;
			zind_sim(t+1)               = markov_draw(zind_sim(t)      ,z_cdf      ,rand_z(t));
		else
			% Draw state tomorrow given innovations
			ssigmaxind_sim(t+1)         = markov_draw(ssigmaxind_sim(t),ssigmax_cdf,rand_unc(t));
			zind_sim(t+1)               = markov_draw(zind_sim(t)      ,z_cdf      ,rand_z(t));
		end
	else
		% Draw state tomorrow given innovations
		ssigmaxind_sim(t+1)             = markov_draw(ssigmaxind_sim(t),ssigmax_cdf,rand_unc(t));
		zind_sim(t+1)                   = markov_draw(zind_sim(t)      ,z_cdf      ,rand_z(t));
	end
end

% Save Results
states_sim.K             = K_sim;
states_sim.dist_k        = dist_k_sim;
states_sim.z             = z_grid(zind_sim)';
states_sim.ssigmax       = ssigmax_grid(ssigmaxind_sim);
states_sim.zind          = zind_sim;
states_sim.ssigmaxind     = ssigmaxind_sim;

controls_sim.C           = C_sim;
controls_sim.w           = w_sim;
controls_sim.q           = q_sim;
controls_sim.ttheta      = ttheta_sim;
controls_sim.inv         = [K_sim(2:T)-(1-ddelta)*K_sim(1:T-1),0];
controls_sim.GDP         = controls_sim.inv.*q_sim + C_sim;
controls_sim.revenue     = revenue_sim;
controls_sim.demand      = demand_sim;
end

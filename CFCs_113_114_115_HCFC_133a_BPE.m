close all
clear 
clc

%% constants
% Samples and years of simulation model (Can be halved for efficiency)
N = 1*10^6;
Nresamps = 0.5*10^6;
Nresamps2 = 2*10^5;
yr_obs_CFC114_CFC115 = 1978:2020;

% molar masses
HCFC133a_mass = 118.5;
HFC134a_mass = 102;
CFC114_mass = 170.9;
CFC113_mass = 187.4;
CFC115_mass = 153.961;

% Conversion factor for emissions/mixing ratios
ppt_to_Gg_CFC113 =  (1.06)/(6.02)*CFC113_mass;
ppt_to_Gg_CFC114 =  (1.06)/(6.02)*CFC114_mass;
ppt_to_Gg_CFC115 =  1.06/6.02*CFC115_mass;
ppt_to_Gg_HCFC133a =  (1.06)/(6.02)*HCFC133a_mass;

% SPARC CFC lifetimes (table 6-1 in Ko et al. 2013; unc. are 2sigma)
lt_CFC113 = 93; % most likely
% lt_CFC113 = 93/0.85; % high end
% lt_CFC113 = 93*0.88; % low end
lt_CFC114 = 191; 
% lt_CFC114 = 191/0.85; 
% lt_CFC114 = 191*0.88;
lt_CFC115 = 540;
% lt_CFC115 = 540/0.85;
% lt_CFC115 = 540*0.88;

% HCFC133a lifetime from Vollmer et al. 2015
lt_HCFC133a = 4.6;

%% HFC production
% HFC-134a production from Velders et al. 2022 (originally from Velders 2015)
% 1990 to 2022, tons per yr (converted to kt/yr)
prod_HFC134a_A5 = 0.001*[0,0,0,0,0,529.2,2097.7,3185.9,4224.8,5031,8502.3, ...
    10636.5,14531.7,19552.2,23873.4,32503.6,38219.2,42947.7,48789,59109.2, ...
    64950.5,77495.9,96597,102365.4,129881,148832.7,159364.4,172779.4,181885.6, ...
    197526.6,214771.4,230396.5,245769.9];
prod_HFC134a_nonA5 = 0.001*[197.2,4398.9,15388.8,31316.5,57818,66192.4,106965.1, ...
    123542.1,138622.7,141050.3,163430.6,151459.9,161943.3,162079.7,160321.3, ...
    160428.2,154975.3,150819.3,158749,174881.1,172247.7,163077.8,163057.5,157795.8, ...
    171290.8,171155.5,160270.3,151751.8,141804.2,112954.7,109027.8,100486.7,96362.8];

% HFC-125 production from Velders 2022 (originally from Velders 2015)
% 1990 to 2022, tons per yr (converted to kt/yr)
prod_HFC125_A5 = 0.001*[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1805.9,2824.1,4352.2, ...
    7331.4,13313.9,35353.6,38883.7,48047.2,58018.1,76917.1,79029.2,92764.1, ...
    105446.2,112007.1,141913.6,157204.8,172071.1];
prod_HFC125_nonA5 = 0.001*[28.7,88.4,62.6,2423.3,6685.5,11769.2,12212.3,16124.7,15837.1,18135, ...
    12375.6,22074.2,31965.1,33287.1,37729.5,47780.8,53009.6,56956.4,63794.6,74069,81487.5,78126.5, ...
    85058.3,76917.8,77106.3,83644.6,82678.4,86177.6,93014.3,75799.1,74875.6,70912.9,70874.6];


%% HCFC-133a emissions and mixing ratios from Vollmer 2021
% Emissions, 1990 to 2022
obs_derived_emiss_HCFC133a = [0.44883,0.35142,0.42542,0.6672,0.59514,0.49215,0.64125,0.50035,0.56548,0.81116,...
    0.33653,0.74621,0.7104,1.02117,1.3762,2.01211,2.13372,2.23916,1.10092,...
    1.48153,2.47457,3.52108,2.53544,1.83928,1.47954,2.3074,2.77171,2.20893,2.15204,1.89945, 1.96790378, 2.315779];

% Mixing ratios, 1990 to 2019
% 1990 to 2019
obs_mixingratios_HCFC133a_archived = [0.04892,0.0573,0.06337,0.07635,0.08943,0.09479,0.10102,0.10543 ...
0.10755,0.11765,0.11752,0.11831,0.12744,0.14175,0.16907,0.21458,0.26578,0.30968,0.31425,0.30454, ...
0.33508,0.40870,0.45876,0.45541,0.43005,0.42787,0.45779,0.47411,0.47240,0.46391];
% 2017 to 2022
obs_mixingratios_HCFC133a_AGAGE = [0.48,0.481,0.466,0.452,0.458,0.4795];
obs_mixingratios_HCFC133a = obs_mixingratios_HCFC133a_archived;
obs_mixingratios_HCFC133a(28:33) = obs_mixingratios_HCFC133a_AGAGE;

%% CFC-115

% from Luke Western's processing of AGAGE data
% 1978 to 2021 obs, 1978 to 2019 emissions
mat = readmatrix("CFC-115_Global_annual_mole_fraction_AGAGE.csv");
CFC115_MF = mat(:,4);
obs_mixingratios_CFC115 = CFC115_MF;

% from Luke Western: either AGAGE or NOAA
mat = readmatrix('CFC-115_Global_annual_emissions_AGAGE.csv');
CFC115_12box_emissions = mat(:,4);
obs_derived_emiss_CFC115 = CFC115_12box_emissions;

% Samples and years of simulation model
simulation_years_CFC115 = 1935:2020;
Nyrs_115 = length(simulation_years_CFC115(1):2020);

% create prior distributions of production for long and short bank uses,
% and release factors from banks
[prod_l_CFC115,prod_s_CFC115, RFlongBank, RFshortBank,Years, short_fraction] = CFC115_priors(N);

% Annual release rate from long banks
rf_l_CFC115 = RFlongBank.LT;
rf_l_CFC115(rf_l_CFC115>0.9) = 0.9;

% Direct emission rate rate during proction for long and short bank uses
de_l_CFC115 = RFlongBank.y1;
de_s_CFC115 = RFshortBank.y1;

% Assume 0 production in years without reported values
prod_l_CFC115(isnan(prod_l_CFC115)) = 0;
prod_s_CFC115(isnan(prod_s_CFC115)) = 0;


% Initiate simulation model
molefractions_y1 = 0; %ppt
t = 1;

% short bank emission: production * direct emission rate
e_s_CFC115(:,t) = prod_s_CFC115(:,t).*de_s_CFC115';
% short bank stroage: production * (1 - direct emission rate)
b_s_CFC115(:,t) = prod_s_CFC115(:,t).*(1 - de_s_CFC115');

% long bank emissions and storage
[b, e] = Bank_Emiss(prod_l_CFC115(:,t), rf_l_CFC115', de_l_CFC115',zeros(N,1));
e_l_CFC115(:,t) = e; 
b_l_CFC115(:,t) = b;

% emissions from HFC-125 production, 1990 to 2019
% byproduction rates for A 5 and non-A 5 countries
byproduct_rate_125_115_A5 = repmat(betarnd(2,2,[1,N])/50,length(prod_HFC125_A5),1)';
byproduct_rate_125_115_nonA5 = repmat(betarnd(2,2,[1,N])/50,length(prod_HFC125_A5),1)'; 
% emission for A 5 and non-A 5 countries (leakage rate * production)
e_HFC125_CFC115_A5 = byproduct_rate_125_115_A5.*prod_HFC125_A5;
e_HFC125_CFC115_nonA5 = byproduct_rate_125_115_nonA5.*prod_HFC125_nonA5;
% global emissions from HFC-125 production
e_HFC125_CFC115 = zeros(N,Nyrs_115);
e_HFC125_CFC115(:,(1990-1934):end) = e_HFC125_CFC115_A5(:,1:31) + e_HFC125_CFC115_nonA5(:,1:31);

% sum emissions
emissions_CFC115(:,t) = e_s_CFC115(:, t) + e_l_CFC115(:, t) + e_HFC125_CFC115(:,t);%  % calculating total emissions from all sources
emissions_CFC115_banks(:,t) = e_s_CFC115(:, t) +  e_l_CFC115(:, t);

% initialize simulation
molefractions_simulation_CFC115_HFC125(:,t+1) = exp(-1./lt_CFC115).*molefractions_y1+(1/ppt_to_Gg_CFC115).*emissions_CFC115(:,t);

% simulate MFs forward with emissions, 1935 to 2020
for t = 2:Nyrs_115 - 1 
    % short banks
    e_s_CFC115(:,t) = prod_s_CFC115(:,t).*de_s_CFC115' + b_s_CFC115(:, t-1);
    b_s_CFC115(:,t) = prod_s_CFC115(:,t).*(1 - de_s_CFC115');
    
    % long banks
    [b, e] = Bank_Emiss(prod_l_CFC115(:,t), rf_l_CFC115', de_l_CFC115',  b_l_CFC115(:,t-1));

    e_l_CFC115(:,t) = e; 
    b_l_CFC115(:,t) = b;
    
    % total emissions
    emissions_CFC115(:,t) = e_s_CFC115(:, t) +  e_l_CFC115(:, t) + e_HFC125_CFC115(:,t); 
    emissions_CFC115_banks(:,t) = e_s_CFC115(:, t) +  e_l_CFC115(:, t);
    
    % "1-box" model (Eq. 1 in methods)
    molefractions_simulation_CFC115_HFC125(:,t+1) = exp(-1./lt_CFC115).*molefractions_simulation_CFC115_HFC125(:,t)+(1/ppt_to_Gg_CFC115)*emissions_CFC115(:,t);
end

% likelihood section 
% years to evaluate likelihood function
year1_likelihood = 1990;
yind1 = find(simulation_years_CFC115 == year1_likelihood);
yind2 = find(simulation_years_CFC115 == 2020);
yind_obs1 = find(yr_obs_CFC114_CFC115 == year1_likelihood); 
yind_obs2 = find(yr_obs_CFC114_CFC115 == 2020); 
n_tmp = length(yind1:yind2);

% uncertainty matrix: 4% of observed values
molefraction_sigma_prior_CFC115 = (0.04*ones(N,1))*obs_mixingratios_CFC115(yind_obs1:yind_obs2)';

% autocorrelation parameter
exp_val = abs(repmat([1:n_tmp],n_tmp,1)-repmat([1:n_tmp]',1,n_tmp));
corrparam = 0.95; % This was selected because we are interested in the long term trend of emissions, not year to year precision
rho_tmp = corrparam*ones(n_tmp,n_tmp);
Rho = rho_tmp.^exp_val;
    
parfor ii = 1:N
    Cov_matrix_CFC115 = Rho.*(molefraction_sigma_prior_CFC115(ii,:)'*molefraction_sigma_prior_CFC115(ii,:));
    likelihood_CFC115_HFC125(ii) = mvnpdf(molefractions_simulation_CFC115_HFC125(ii, yind1:yind2), obs_mixingratios_CFC115(yind_obs1:yind_obs2)', Cov_matrix_CFC115);
end

Likelihood_CFC115_HFC125 = (1/sum(likelihood_CFC115_HFC125))*likelihood_CFC115_HFC125;

figure; 
plot((1/sum(Likelihood_CFC115_HFC125))*cumsum(Likelihood_CFC115_HFC125))
title('Important Ratio CFC-115');


% Resample from the priors based on the relative likelihood
ResampleIndex_CFC115= nan(Nresamps,1);

% Updating priors based on sampling importantce ratio/likelihood 
IR_vec_HFC125 = 1/sum(Likelihood_CFC115_HFC125)*cumsum(Likelihood_CFC115_HFC125);

% identifying samples that will be resampled in the posterior
for ii = 1:Nresamps
    
    ResampleIndex_CFC115(ii) = find(IR_vec_HFC125>rand,1);
end 


%% CFC-114
% steps are the same as those described in the CFC-115 section
N = Nresamps;

% Mixing ratios from Western's output
mat = readmatrix("CFC-114_Global_annual_mole_fraction_AGAGE.csv");
obs_mixingratios_CFC114 = mat(:,4);

% emissions from Western's output
mat = readmatrix('CFC-114_Global_annual_emissions_AGAGE.csv');
obs_derived_emiss_CFC114 = mat(:,4);


% Initiate simulation 
simulation_years = 1935:2020;
simulation_years_CFC114 = simulation_years;
Nyrs_114 = length(simulation_years(1):2020);
t = 1;

[prod_l_CFC114, prod_s_CFC114, RFlongBank, RFshortBank, Years] = CFC114_priors(N);
rf_l_CFC114 = RFlongBank.LT;
de_l_CFC114 = RFlongBank.y1; 
de_s_CFC114 = RFshortBank.y1; 

t = 1;
[b, e] = Bank_Emiss(prod_l_CFC114(:,t), rf_l_CFC114', de_l_CFC114',zeros(N,1));
e_l_CFC114(:,t) = e; 
b_l_CFC114(:,t) = b;

e_s_CFC114(:,t) = prod_s_CFC114(:,t).*de_s_CFC114';
b_s_CFC114(:,t) = prod_s_CFC114(:,t).*(1 - de_s_CFC114');

% chemical conversion rate of 114 -> 134a
conversion_rate_114_134a_nonA5 = 0.94; %  0.94 from Morikawa 1994 patent
conversion_rate_114_134a_A5 = conversion_rate_114_134a_nonA5;

% factor allocating between 114 pathway and 133a pathway
% Pathway is time varying (can go both up and down, but 
% with high correlation so the changes aren't too abrupt)
% A 5 countries
Nyrs_133a = 33; 
rho(1,1,:) = 0.95+0.05*betarnd(2,2,N,1);
rho_tmp = repmat(rho,Nyrs_133a,Nyrs_133a,1);
exp_val = repmat(abs(repmat([1:Nyrs_133a],Nyrs_133a,1)-repmat([1:Nyrs_133a]',1,Nyrs_133a)),1,1,N);
Rhom = rho_tmp.^exp_val;
clear rho_tmp exp_val
Zm = mvnrnd(zeros(N,Nyrs_133a), Rhom, N);
clear Rhom
Um = normcdf(Zm,0,1);

% At most 70% of the production in A5 is coming from CFC-114 pathway
pathway_fraction_A5 = 0.7*Um;

% Non-A 5 countries
Nyrs_133a = 33; 
rho(1,1,:) = 0.95+0.05*betarnd(2,2,N,1);
rho_tmp = repmat(rho,Nyrs_133a,Nyrs_133a,1);
exp_val = repmat(abs(repmat([1:Nyrs_133a],Nyrs_133a,1)-repmat([1:Nyrs_133a]',1,Nyrs_133a)),1,1,N);
Rhom = rho_tmp.^exp_val;
clear rho_tmp exp_val
Zm = mvnrnd(zeros(N,Nyrs_133a), Rhom, N);
clear Rhom
Um = normcdf(Zm,0,1);
% At most 70% of the production in non-A5 is coming from CFC-114 pathway
pathway_fraction_nonA5 = 0.7*Um;

% A 5 CFC-114 production based on HFC-134a production 
emission_rate_CFC114_FS_A5 = repmat(betarnd(2,2,[1,N])/50,33,1)';
A5_114_FS_prod = pathway_fraction_A5.*prod_HFC134a_A5(1:33)./(conversion_rate_114_134a_A5*(1-emission_rate_CFC114_FS_A5))*(CFC114_mass/HFC134a_mass);
% non-A 5 
emission_rate_CFC114_FS_nonA5 = repmat(betarnd(2,2,[1,N])/50,33,1)';
nonA5_114_FS_prod = pathway_fraction_nonA5.*prod_HFC134a_nonA5(1:33)./(conversion_rate_114_134a_nonA5*(1-emission_rate_CFC114_FS_nonA5))*(CFC114_mass/HFC134a_mass);

% A5 CFC-114 emissions from feedstock production for HFC-134a 
% feedstock emission rate: beta dist btwn 0 and 2%
e_nonA5_FS_CFC114 = zeros(N,88);
e_nonA5_FS_CFC114(:,(1990-1934):(2022-1934)) = nonA5_114_FS_prod.*emission_rate_CFC114_FS_nonA5;
e_nonA5_FS_CFC114(e_nonA5_FS_CFC114 < 0) = 0;

% A5 CFC-114 emissions from feedstock production for HFC-134a 
emission_rate_CFC114_FS_A5 = repmat(betarnd(2,2,[1,N])/50,33,1)';
e_A5_FS_CFC114 = zeros(N,88);
e_A5_FS_CFC114(:,(1990-1934):(2022-1934)) = A5_114_FS_prod.*emission_rate_CFC114_FS_A5;
e_A5_FS_CFC114(e_A5_FS_CFC114 < 0) = 0;

% emissions from 125 production, 1990 to 2020
% Leakage rate for A 5 and non-A 5 
CFC114to115_byprod = 0.5*rand(N,1);
byproduct_rate_125_114_A5 = CFC114to115_byprod.*byproduct_rate_125_115_A5(ResampleIndex_CFC115,:); 
byproduct_rate_125_114_nonA5 = CFC114to115_byprod.*byproduct_rate_125_115_nonA5(ResampleIndex_CFC115,:);
% emissions: byprod rate * production
e_HFC125_CFC114_A5 = byproduct_rate_125_114_A5.*prod_HFC125_A5;
e_HFC125_CFC114_nonA5 = byproduct_rate_125_114_nonA5.*prod_HFC125_nonA5;
% global emissions from HFC-125 production (A 5 + non-A 5)
e_HFC125_CFC114 = zeros(N,Nyrs_114);
e_HFC125_CFC114(:,(1990-1934):end) = e_HFC125_CFC114_A5(:,1:31) + e_HFC125_CFC114_nonA5(:,1:31);


%% HCFC-133a, which depends on CFC-114 usasge in HFC-134a production
% divide HFC-134a production by pathway 
pathway_1_HFC134a_nonA5 = pathway_fraction_nonA5.*repmat(prod_HFC134a_nonA5, N, 1);
pathway_2_HFC134a_nonA5 = (1-pathway_fraction_nonA5).*repmat(prod_HFC134a_nonA5, N, 1);
pathway_1_HFC134a_A5 = pathway_fraction_A5.*repmat(prod_HFC134a_A5,N,1);
pathway_2_HFC134a_A5 = (1-pathway_fraction_A5).*repmat(prod_HFC134a_A5, N, 1);

conversion_rate_133a_134a = 0.95; % From Scott and Steve

% emission rate of HCFC-133a in non-A 5 and A 5 countries from HFC-134a
% beta dist btwn 0 and 1/66
emission_rate_HFC134a_HCFC133a_nonA5 = repmat(betarnd(2,2,[1,N])/66,33,1)'; 
emission_rate_HFC134a_HCFC133a_A5 = repmat(betarnd(2,2,[1,N])/66,33,1)'; 


% HCFC133a production = HFC134a prod/Conversion rate
HCFC133a_HFC134a_nonA5 = pathway_2_HFC134a_nonA5./(conversion_rate_133a_134a*(1-emission_rate_HFC134a_HCFC133a_nonA5))*(HCFC133a_mass/HFC134a_mass);
HCFC133a_HFC134a_A5 = pathway_2_HFC134a_A5./(conversion_rate_133a_134a*(1-emission_rate_HFC134a_HCFC133a_A5))*(HCFC133a_mass/HFC134a_mass);

% emissions of HCFC-133a in non-A 5 and A 5 countries
e_HFC134a_HCFC133a_nonA5 = HCFC133a_HFC134a_nonA5.*emission_rate_HFC134a_HCFC133a_nonA5;
e_HFC134a_HCFC133a_A5 = HCFC133a_HFC134a_A5.*emission_rate_HFC134a_HCFC133a_A5;

% emission rate of HCFC-133a in non-A 5 and A 5 countries from HFC-125
% beta dist btwn 0 and 1/66
emission_rate_HFC125_HCFC133a_nonA5 = repmat(betarnd(2,2,[1,N])/66,33,1)'; 
emission_rate_HFC125_HCFC133a_A5 = repmat(betarnd(2,2,[1,N])/66,33,1)'; 

% emissions of HCFC-133a in non-A 5 and A 5 countries
e_HFC125_HCFC133a_nonA5 = prod_HFC125_nonA5(1:33).*emission_rate_HFC125_HCFC133a_nonA5;
e_HFC125_HCFC133a_A5 = prod_HFC125_A5(1:33).*emission_rate_HFC125_HCFC133a_A5;

%% simulate mixing ratios of each compound (same steps as CFC115 section)
% Initiating simulation model
molefractions_y1 = 0; %ppt
t = 1; 

% emissions for first time step
emissions_CFC114_total(:,t) = e_s_CFC114(:, t) + e_l_CFC114(:, t) + e_nonA5_FS_CFC114(:, t) + e_A5_FS_CFC114(:,t) + e_HFC125_CFC114(:,t); 
emissions_HCFC133a_total = e_HFC134a_HCFC133a_nonA5 + e_HFC134a_HCFC133a_A5 + e_HFC125_HCFC133a_nonA5 + e_HFC125_HCFC133a_A5;

% initialize mole fractions
molefractions_simulation_CFC114_total(:,t+1) = exp(-1./lt_CFC114).*molefractions_y1+(1/ppt_to_Gg_CFC114).*emissions_CFC114_total(:,t);

molefractions_HCFC133a_y1 = 0.0489; % mixing ratio in 1990
molefractions_simulation_HCFC133a(:,1) = molefractions_HCFC133a_y1*ones(N,1);
molefractions_simulation_HCFC133a(:,t+1) =  exp(-1./lt_HCFC133a).*molefractions_HCFC133a_y1+(1/ppt_to_Gg_HCFC133a).*emissions_HCFC133a_total(:,t);

% simulate mole fractions, CFC-114
for t = 2:Nyrs_114 - 1
    e_s_CFC114(:,t) = prod_s_CFC114(:,t).*de_s_CFC114' + b_s_CFC114(:, t-1);
    b_s_CFC114(:,t) = prod_s_CFC114(:,t).*(1 - de_s_CFC114');

    [b, e] = Bank_Emiss(prod_l_CFC114(:,t), rf_l_CFC114', de_l_CFC114',  b_l_CFC114(:,t-1));

    e_l_CFC114(:,t) = e; 
    b_l_CFC114(:,t) = b;

    emissions_CFC114_total(:,t) = e_s_CFC114(:, t) + e_l_CFC114(:, t) + e_nonA5_FS_CFC114(:, t) + e_A5_FS_CFC114(:,t) + e_HFC125_CFC114(:,t);

    molefractions_simulation_CFC114_total(:,t+1) = exp(-1./lt_CFC114).*molefractions_simulation_CFC114_total(:,t)+(1/ppt_to_Gg_CFC114)*emissions_CFC114_total(:,t);

end

% simulate mole fractions, HCFC-133a
for t = 2:(Nyrs_133a)-1
    molefractions_simulation_HCFC133a(:,t+1) =  exp(-1./lt_HCFC133a).*molefractions_simulation_HCFC133a(:,t)+(1/ppt_to_Gg_HCFC133a).*emissions_HCFC133a_total(:,t);
end

%% likelihood section
year1_likelihood = 1990;

yind1_114 = find(simulation_years_CFC114 == year1_likelihood);
yind2_114 = find(simulation_years_CFC114 == 2020);

yind_obs1 = find(yr_obs_CFC114_CFC115 == year1_likelihood); 
yind_obs2 = find(yr_obs_CFC114_CFC115 == 2020); 

n_tmp = length(yind1_114:yind2_114);

abs_uncertainty_133a = 0.2*ones(N,1); 
molefraction_sigma_prior_HCFC133a = abs_uncertainty_133a*ones(1,6); 
Obs_5yrs = mean(reshape(obs_mixingratios_HCFC133a(1:30),5,6),1);

molefraction_sigma_prior_CFC114_uncertainty = 0.04*ones(N,length(yind_obs1:yind_obs2));
molefraction_sigma_prior_CFC114 = (molefraction_sigma_prior_CFC114_uncertainty).*obs_mixingratios_CFC114(yind_obs1:yind_obs2)';

exp_val = abs(repmat([1:n_tmp],n_tmp,1)-repmat([1:n_tmp]',1,n_tmp));
 
corrparam_114 = 0.95; 
corrparam_133a = 0.6 + 0.2*betarnd(2,2,N,1); 

rho_tmp = corrparam_114*ones(n_tmp,n_tmp);
Rho_114 =  rho_tmp.^exp_val;

n134a = size(molefraction_sigma_prior_HCFC133a, 2);
exp_val = abs(repmat([1:n134a],n134a,1)-repmat([1:n134a]',1,n134a));
rho_tmp = corrparam_114*ones(n134a,n134a);
Rho_133a = rho_tmp.^exp_val;

for ii = 1:N
    Cov_matrix_CFC114 = Rho_114.*(molefraction_sigma_prior_CFC114(ii,:)'*molefraction_sigma_prior_CFC114(ii,:));
    Cov_matrix_HCFC133a = Rho_133a.*(molefraction_sigma_prior_HCFC133a(ii,:)'*molefraction_sigma_prior_HCFC133a(ii,:));
    
    likelihood_CFC114_total(ii) = mvnpdf(molefractions_simulation_CFC114_total(ii, yind1_114:yind2_114), obs_mixingratios_CFC114(yind_obs1:yind_obs2)', Cov_matrix_CFC114);
    likelihood_HCFC133a(ii) = mvnpdf(molefractions_simulation_HCFC133a(ii, 3:5:30), Obs_5yrs, Cov_matrix_HCFC133a);

end

Likelihood_CFC114 = (1/sum(likelihood_CFC114_total))*cumsum(likelihood_CFC114_total);
Likelihood_HCFC133a = (1/sum(likelihood_HCFC133a))*cumsum(likelihood_HCFC133a);

likelihood_CFC114_HCFC133a_joint = likelihood_CFC114_total.*likelihood_HCFC133a;
Likelihood_CFC114_HCFC133a_Joint = (1/sum(likelihood_CFC114_HCFC133a_joint))*cumsum(likelihood_CFC114_HCFC133a_joint);

% Checking the importance ratio is relatively linear
figure
subplot(1,3,1)
plot(Likelihood_CFC114)
subplot(1,3,2)
plot(Likelihood_HCFC133a)
subplot(1,3,3)
plot(Likelihood_CFC114_HCFC133a_Joint)


%% Update priors based on sampling importantce ratio/likelihood 

IR_vec_joint = Likelihood_CFC114_HCFC133a_Joint;

ResampleIndex_CFC114_HCFC133a_Joint = nan(Nresamps,1);

for ii = 1:Nresamps
    ResampleIndex_CFC114_HCFC133a_Joint(ii) = find(IR_vec_joint>rand,1);
end 

%% Sequential updatig from 114/133a IR

indx_resample = ResampleIndex_CFC114_HCFC133a_Joint;


%% CFC-113, including conversion from CFC-114
% same steps as previous sections

% Mixing ratios, from Western again
mat = readmatrix("CFC-113_Global_annual_mole_fraction_AGAGE.csv");
years_12box = mat(:,1);
obs_mixingratios_CFC113 = mat(:,4);

% observed emissions, from Western
mat = readmatrix('CFC-113_Global_annual_emissions_AGAGE.csv');
obs_derived_emiss_CFC113 = mat(:,4);

% Samples and years of simulation model
simulation_years_CFC113 = 1955:2020;
Nyrs_113 = length(simulation_years_CFC113);

[prod_l_CFC113,prod_s_CFC113,RFlongBankLT, de_l_CFC113, de_s_CFC113,Years] = CFC113_priors(Nresamps);
prod_l_CFC113 = 0.001*prod_l_CFC113;
prod_s_CFC113 = 0.001*prod_s_CFC113;

rf_l_CFC113 = 1 - exp(-1./RFlongBankLT);

prod_s_CFC113(:,65:68) = 0;


% emissions from 125 production, 1990 to 2022
% Byprod rate for  A 5 and non-A 5
CFC113to115_byprod = 0.5*rand(Nresamps,1);
byproduct_rate_125_113_A5 = CFC113to115_byprod.*byproduct_rate_125_115_A5(ResampleIndex_CFC115,:); 
byproduct_rate_125_113_nonA5 = CFC113to115_byprod.*byproduct_rate_125_115_nonA5(ResampleIndex_CFC115,:);
% emission (leakage rate * production)
e_HFC125_CFC113_A5 = byproduct_rate_125_113_A5.*prod_HFC125_A5;
e_HFC125_CFC113_nonA5 = byproduct_rate_125_113_nonA5.*prod_HFC125_nonA5;
% global emissions from HFC-125 production
e_HFC125_CFC113 = zeros(Nresamps,Nyrs_113);
e_HFC125_CFC113(:,(1990-1954):end) = e_HFC125_CFC113_A5(:,1:31) + e_HFC125_CFC113_nonA5(:,1:31);


%%%%%%%%%%%%%%% THIS DEPENDS ON CFC-114 %%%%%%%%%%%%%%%%%%%
% non-A 5 and A 5 CFC-114 and CFC-113 production
nonA5_114_FS_prod_constrained = nonA5_114_FS_prod(indx_resample,:);
A5_114_FS_prod_constrained = A5_114_FS_prod(indx_resample,:);

% 113 to 114a conversion rate
conversion_rate_113_114_nonA5 = 0.98; % from Grumprecht 1991 patent
conversion_rate_113_114_A5 = conversion_rate_113_114_nonA5; 

% non-A 5 CFC-113 feedstock production
emission_rate_CFC113_FS_nonA5 = repmat(betarnd(2,2,[1,Nresamps])/25,33,1)';
nonA5_113_FS_prod = (1-emission_rate_CFC113_FS_nonA5).*nonA5_114_FS_prod_constrained./conversion_rate_113_114_nonA5*(CFC113_mass/CFC114_mass);
% A 5 CFC-113 feedstock production
emission_rate_CFC113_FS_A5 = repmat(betarnd(2,2,[1,Nresamps])/25,33,1)';
A5_113_FS_prod = (1-emission_rate_CFC113_FS_A5).*A5_114_FS_prod_constrained./conversion_rate_113_114_A5*(CFC113_mass/CFC114_mass);

total_113_FS_prod  = nonA5_113_FS_prod + A5_113_FS_prod;
total_114_FS_prod  = nonA5_114_FS_prod_constrained + A5_114_FS_prod_constrained;

% non-A 5 CFC-113 feedstock emissions
e_nonA5_FS_CFC113 = zeros(Nresamps,68);
e_nonA5_FS_CFC113(:,(1990-1954):(2022-1954)) = nonA5_113_FS_prod.*emission_rate_CFC113_FS_nonA5;
e_nonA5_FS_CFC113(e_nonA5_FS_CFC113 < 0) = 0;

% A 5 CFC-113 feedstock emissions
e_A5_FS_CFC113 = zeros(Nresamps,68);
e_A5_FS_CFC113(:,(1990-1954):(2022-1954)) = A5_113_FS_prod.*emission_rate_CFC113_FS_A5;
e_A5_FS_CFC113(e_A5_FS_CFC113 < 0) = 0;


% Initiate simulation model
molefractions_y1 = 0; %ppt
t = 1;

e_s_CFC113(:,t) = prod_s_CFC113(:,t).*de_s_CFC113';
b_s_CFC113(:,t) = prod_s_CFC113(:,t).*(1 - de_s_CFC113');

prod_l_CFC113(:,65:68) = 0;

[b, e] = Bank_Emiss(prod_l_CFC113(:,t), rf_l_CFC113', de_l_CFC113',zeros(Nresamps,1));
e_l_CFC113(:,t) = e; 
b_l_CFC113(:,t) = b;

% % emissions for first time step
emissions_CFC113(:,t) = e_s_CFC113(:, t) + e_l_CFC113(:, t) + e_nonA5_FS_CFC113(:, t) + e_A5_FS_CFC113(:,t) + e_HFC125_CFC113(:,t); % emissions from all sources 
molefractions_simulation_CFC113(:,t+1) = exp(-1./lt_CFC113).*molefractions_y1+(1/ppt_to_Gg_CFC113).*emissions_CFC113(:,t);

% simulate mole fractions, CFC-113
for t = 2:Nyrs_113 - 1
    e_s_CFC113(:,t) = prod_s_CFC113(:,t).*de_s_CFC113' + b_s_CFC113(:, t-1);
    b_s_CFC113(:,t) = prod_s_CFC113(:,t).*(1 - de_s_CFC113');

    [b, e] = Bank_Emiss(prod_l_CFC113(:,t), rf_l_CFC113', de_l_CFC113',  b_l_CFC113(:,t-1));

    e_l_CFC113(:,t) = e; 
    b_l_CFC113(:,t) = b;

    emissions_CFC113(:,t) = e_s_CFC113(:, t) + e_l_CFC113(:, t) + e_nonA5_FS_CFC113(:, t) + e_A5_FS_CFC113(:,t)  + e_HFC125_CFC113(:,t);
    molefractions_simulation_CFC113(:,t+1) = exp(-1./lt_CFC113).*molefractions_simulation_CFC113(:,t)+(1/ppt_to_Gg_CFC113)*emissions_CFC113(:,t);
end


%% likelihood section
year1_likelihood = 1990;
yind1_113 = find(simulation_years_CFC113 == year1_likelihood);
yind2_113 = find(simulation_years_CFC113 == 2020);
yr_obs_CFC113 = 1982:2020;

yind_obs1 = find(yr_obs_CFC113 == year1_likelihood); 
yind_obs2 = find(yr_obs_CFC113 == 2020); 

molefraction_sigma_prior_CFC113_uncertainty = 0.03*ones(Nresamps,length(yind_obs1:yind_obs2));
molefraction_sigma_prior_CFC113 = (molefraction_sigma_prior_CFC113_uncertainty).*obs_mixingratios_CFC113(yind_obs1:yind_obs2)';

exp_val = abs(repmat([1:n_tmp],n_tmp,1)-repmat([1:n_tmp]',1,n_tmp));

corrparam_113 = 0.95; 

rho_tmp = corrparam_113*ones(n_tmp,n_tmp);
Rho_113 = rho_tmp.^exp_val;


for ii = 1:Nresamps
    Cov_matrix_CFC113 = Rho_113.*(molefraction_sigma_prior_CFC113(ii,:)'*molefraction_sigma_prior_CFC113(ii,:));
    likelihood_CFC113_nonA5_A5(ii) = mvnpdf(molefractions_simulation_CFC113(ii, yind1_113:yind2_113), obs_mixingratios_CFC113(yind_obs1:yind_obs2)', Cov_matrix_CFC113);
end


Likelihood_CFC113 = (1/sum(likelihood_CFC113_nonA5_A5))*cumsum(likelihood_CFC113_nonA5_A5);

figure 
plot(Likelihood_CFC113)

% Update priors based on sampling importantce ratio/likelihood 
IR_vec_CFC113 = Likelihood_CFC113;

ResampleIndex_CFC113 = nan(Nresamps2,1);

parfor ii = 1:Nresamps2
    ResampleIndex_CFC113(ii) = find(IR_vec_CFC113>rand,1);
end 

ResampleIndex_all_joint = indx_resample(ResampleIndex_CFC113);

%% calculate mixing ratios from 2004 to 2020 with banks only
molefractions_simulation_CFC113_banks = zeros(length(ResampleIndex_CFC113),17);
molefractions_simulation_CFC114_banks = zeros(length(ResampleIndex_all_joint),17);
molefractions_simulation_CFC115_banks = zeros(length(ResampleIndex_CFC115),17);
molefractions_simulation_HCFC133a_banks = zeros(length(ResampleIndex_all_joint),17);

molefractions_simulation_CFC113_banks(:,1) = molefractions_simulation_CFC113(ResampleIndex_CFC113,end-16); 
molefractions_simulation_CFC114_banks(:,1) = molefractions_simulation_CFC114_total(ResampleIndex_all_joint,end-16); 
molefractions_simulation_CFC115_banks(:,1) = molefractions_simulation_CFC115_HFC125(ResampleIndex_CFC115,end-16); 
molefractions_simulation_HCFC133a_banks(:,1) = molefractions_simulation_HCFC133a(ResampleIndex_all_joint,end-18); 

e_banks_CFC113 = e_l_CFC113 + e_s_CFC113;
e_banks_CFC114 = e_l_CFC114 + e_s_CFC114;
e_banks_CFC115 = e_l_CFC115 + e_s_CFC115;
for t = 2:17
    molefractions_simulation_CFC113_banks(:,t) = exp(-1./lt_CFC113).*molefractions_simulation_CFC113_banks(:,t-1)+(1/ppt_to_Gg_CFC113)*(e_banks_CFC113(ResampleIndex_CFC113, end-17+t));
    molefractions_simulation_CFC114_banks(:,t) = exp(-1./lt_CFC114).*molefractions_simulation_CFC114_banks(:,t-1)+(1/ppt_to_Gg_CFC114)*(e_banks_CFC114(ResampleIndex_all_joint, end-17+t));
    molefractions_simulation_CFC115_banks(:,t) = exp(-1./lt_CFC115).*molefractions_simulation_CFC115_banks(:,t-1)+(1/ppt_to_Gg_CFC115)*(e_banks_CFC115(ResampleIndex_CFC115, end-17+t));
    molefractions_simulation_HCFC133a_banks(:,t) = exp(-1./lt_HCFC133a).*molefractions_simulation_HCFC133a_banks(:,t-1);
end
% 
% plot(mean(molefractions_simulation_CFC113_banks))
% hold on
% plot(mean(molefractions_simulation_CFC114_banks))
% plot(mean(molefractions_simulation_CFC115_banks))


%% Mixing ratios
fig = figure;

% Figure 3
subplot(4,2,1)
LB = prctile(molefractions_simulation_CFC113_banks,16);
MED = prctile(molefractions_simulation_CFC113_banks,50);
UB  = prctile(molefractions_simulation_CFC113_banks,84);
p2 = boundedline(2004:2020, MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0. 0. 0.7]);
hold on
LB = prctile(molefractions_simulation_CFC113(ResampleIndex_CFC113,:),16);
MED = prctile(molefractions_simulation_CFC113(ResampleIndex_CFC113,:),50);
UB  = prctile(molefractions_simulation_CFC113(ResampleIndex_CFC113,:),84);
p3 = boundedline(1955:2020, MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0.7 0 0]);
p1 = plot(1982:2021,obs_mixingratios_CFC113,'k','LineWidth',2,'LineStyle','--');
xlim([2004,2019])
xticks([])
ylabel({'CFC-113','Mixing ratios (ppt)'})
leg = legend([p1,p2,p3],'Obs.','No HFC production','Posterior','NumColumns',2);
leg.ItemTokenSize(1) = 18;
legend boxoff 
set(gca,'FontSize',12)
set(legend,'FontSize',12)
set(gca,'FontSize',12)

% CFC-114 mixing ratios
subplot(4,2,3)
LB = prctile(molefractions_simulation_CFC114_banks,16);
MED = prctile(molefractions_simulation_CFC114_banks,50);
UB  = prctile(molefractions_simulation_CFC114_banks,84);
p2 = boundedline(2004:2020, MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0. 0. 0.7]);
hold on
LB = prctile(molefractions_simulation_CFC114_total(ResampleIndex_all_joint,:),16);
MED = prctile(molefractions_simulation_CFC114_total(ResampleIndex_all_joint,:),50);
UB  = prctile(molefractions_simulation_CFC114_total(ResampleIndex_all_joint,:),84);
p3 = boundedline(1935:2020, MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0.7 0 0]);
plot(1978:2021,obs_mixingratios_CFC114,'k','LineWidth',2,'LineStyle','--')
xlim([2004,2019])
xticks([])
ylim([15.2,17.2])
% title('CFC-114')
ylabel({'CFC-114','Mixing ratios (ppt)'})
% ylabel('Mixing ratios (ppt)')
set(gca,'FontSize',12)

% CFC-115 mixing ratios
subplot(4,2,5)
LB = prctile(molefractions_simulation_CFC115_banks,16);
MED = prctile(molefractions_simulation_CFC115_banks,50);
UB  = prctile(molefractions_simulation_CFC115_banks,84);
p2 = boundedline(2004:2020, MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0. 0. 0.7]);
hold on
LB = prctile(molefractions_simulation_CFC115_HFC125(ResampleIndex_CFC115,:),16);
MED = prctile(molefractions_simulation_CFC115_HFC125(ResampleIndex_CFC115,:),50);
UB  = prctile(molefractions_simulation_CFC115_HFC125(ResampleIndex_CFC115,:),84);
p3 = boundedline(1935:2020, MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0.7 0 0]);
plot(1978:2021,obs_mixingratios_CFC115,'k','LineWidth',2,'LineStyle','--')
xlim([2004,2019])
xticks([])
ylabel({'CFC-115','Mixing ratios (ppt)'})
set(gca,'FontSize',12)

% HCFC-133a mixing ratios
subplot(4,2,7)
LB = prctile(molefractions_simulation_HCFC133a_banks,16);
MED = prctile(molefractions_simulation_HCFC133a_banks,50);
UB  = prctile(molefractions_simulation_HCFC133a_banks,84);
p2 = boundedline(2004:2020, MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0. 0. 0.7]);
hold on
LB = prctile(molefractions_simulation_HCFC133a(ResampleIndex_all_joint,:),16);
MED = prctile(molefractions_simulation_HCFC133a(ResampleIndex_all_joint,:),50);
UB  = prctile(molefractions_simulation_HCFC133a(ResampleIndex_all_joint,:),84);
p3 = boundedline(1990:2022, MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0.7 0 0]);
plot(1990:2022,obs_mixingratios_HCFC133a,'k','LineWidth',2,'LineStyle','--')
xlim([2004,2019])
ylabel({'HCFC-133a','Mixing ratios (ppt)'})
set(gca,'FontSize',12)
xlabel('Year')

% Emissions and Emissions by source
% CFC-113
subplot(4,2,2)
e_FS_CFC113 = e_A5_FS_CFC113 + e_nonA5_FS_CFC113;
hold on
LB = prctile(emissions_CFC113(ResampleIndex_CFC113,:),16);
MED = prctile(emissions_CFC113(ResampleIndex_CFC113,:),50);
UB  = prctile(emissions_CFC113(ResampleIndex_CFC113,:),84);
p1 = boundedline(1955:2019, MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0.7 0 0]);
LB = prctile(e_banks_CFC113(ResampleIndex_CFC113,:),16);
MED = prctile(e_banks_CFC113(ResampleIndex_CFC113,:),50);
UB  = prctile(e_banks_CFC113(ResampleIndex_CFC113,:),84);
p2 = boundedline(1955:2019, MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0 0. 0.7]);
hold on
LB = prctile(e_FS_CFC113(ResampleIndex_CFC113,:),16);
MED = prctile(e_FS_CFC113(ResampleIndex_CFC113,:),50);
UB  = prctile(e_FS_CFC113(ResampleIndex_CFC113,:),84);
p3 = boundedline(1955:2022, MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0 0.7 0.]);
LB = prctile(e_HFC125_CFC113(ResampleIndex_CFC113,:),16);
MED = prctile(e_HFC125_CFC113(ResampleIndex_CFC113,:),50);
UB  = prctile(e_HFC125_CFC113(ResampleIndex_CFC113,:),84);
p4 = boundedline(1955:2020, MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0.7 0 0.7]);
p0 = plot(1982:2021,obs_derived_emiss_CFC113,'k','LineWidth',2,'LineStyle','--');
xlim([2004,2019])
xticks([])
ylim([0,15])
ylabel('Emissions (Gg\cdoty^{-1})')
leg = legend([p0,p1,p2,p3,p4], 'Obs. derived','Posterior total','Production/Banks','HFC-134a','HFC-125','NumColumns',2);
leg.ItemTokenSize(1) = 18;
legend boxoff 
set(gca,'FontSize',12)
set(legend,'FontSize',12)

% CFC-114
subplot(4,2,4)
e_FS_CFC114 = e_A5_FS_CFC114 + e_nonA5_FS_CFC114;
hold on
LB = prctile(emissions_CFC114_total(ResampleIndex_all_joint,:),16);
MED = prctile(emissions_CFC114_total(ResampleIndex_all_joint,:),50);
UB  = prctile(emissions_CFC114_total(ResampleIndex_all_joint,:),84);
p1 = boundedline(1935:2019, MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0.7 0 0]);
LB = prctile(e_banks_CFC114(ResampleIndex_all_joint,:),16);
MED = prctile(e_banks_CFC114(ResampleIndex_all_joint,:),50);
UB  = prctile(e_banks_CFC114(ResampleIndex_all_joint,:),84);
p2 = boundedline(1935:2019, MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0. 0. 0.7]);
hold on
LB = prctile(e_FS_CFC114(ResampleIndex_all_joint,:),16);
MED = prctile(e_FS_CFC114(ResampleIndex_all_joint,:),50);
UB  = prctile(e_FS_CFC114(ResampleIndex_all_joint,:),84);
p3 = boundedline(1935:2022, MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0 0.7 0.]);
LB = prctile(e_HFC125_CFC114(ResampleIndex_all_joint,:),16);
MED = prctile(e_HFC125_CFC114(ResampleIndex_all_joint,:),50);
UB  = prctile(e_HFC125_CFC114(ResampleIndex_all_joint,:),84);
p4 = boundedline(1935:2020, MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0.7 0 0.7]);
p0 = plot(1978:2021,obs_derived_emiss_CFC114,'k','LineWidth',2,'LineStyle','--');
xlim([2004,2019])
xticks([])
ylim([0,4])
ylabel('Emissions (Gg\cdoty^{-1})')
set(gca,'FontSize',12)

% CFC-115
subplot(4,2,6)
LB = prctile(emissions_CFC115(ResampleIndex_CFC115,:),16);
MED = prctile(emissions_CFC115(ResampleIndex_CFC115,:),50);
UB  = prctile(emissions_CFC115(ResampleIndex_CFC115,:),84);
p2 = boundedline(1935:2019, MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0.7 0 0]);
hold on
LB = prctile(e_banks_CFC115(ResampleIndex_CFC115,:),16);
MED = prctile(e_banks_CFC115(ResampleIndex_CFC115,:),50);
UB  = prctile(e_banks_CFC115(ResampleIndex_CFC115,:),84);
p1 = boundedline(1935:2019, MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0. 0 0.7]);
LB = prctile(e_HFC125_CFC115(ResampleIndex_CFC115,:),16);
MED = prctile(e_HFC125_CFC115(ResampleIndex_CFC115,:),50);
UB  = prctile(e_HFC125_CFC115   (ResampleIndex_CFC115,:),84);
p2 = boundedline(1935:2020, MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0.7 0 0.7]);
p0 = plot(1978:2021,obs_derived_emiss_CFC115,'k','LineWidth',2,'LineStyle','--');
xlim([2004,2019])
xticks([])
ylabel('Emissions (Gg\cdoty^{-1})')
set(gca,'FontSize',12)

% HCFC-133a
emissions_HCFC133a_HFC134a = e_HFC134a_HCFC133a_A5 + e_HFC134a_HCFC133a_nonA5;
emissions_HCFC133a_HFC125 = e_HFC125_HCFC133a_A5 + e_HFC125_HCFC133a_nonA5;
subplot(4,2,8)
LB = prctile(emissions_HCFC133a_total(ResampleIndex_all_joint,:),16);
MED = prctile(emissions_HCFC133a_total(ResampleIndex_all_joint,:),50);
UB  = prctile(emissions_HCFC133a_total(ResampleIndex_all_joint,:),84);
p2 = boundedline(1990:2022, MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0.7 0 0]);
hold on
LB = prctile(emissions_HCFC133a_HFC134a(ResampleIndex_all_joint,:),16);
MED = prctile(emissions_HCFC133a_HFC134a(ResampleIndex_all_joint,:),50);
UB  = prctile(emissions_HCFC133a_HFC134a(ResampleIndex_all_joint,:),84);
p3 = boundedline(1990:2022, MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0. 0.7 0]);
LB = prctile(emissions_HCFC133a_HFC125(ResampleIndex_all_joint,:),16);
MED = prctile(emissions_HCFC133a_HFC125(ResampleIndex_all_joint,:),50);
UB  = prctile(emissions_HCFC133a_HFC125(ResampleIndex_all_joint,:),84);
p3 = boundedline(1990:2022, MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0.7 0. 0.7]);
plot(1990:2021,obs_derived_emiss_HCFC133a,'k','LineWidth',2,'LineStyle','--')
xlim([2004,2019])
ylabel('Emissions (Gg\cdoty^{-1})')
set(gca,'FontSize',12)
xlabel('Year')

figure_width = 8; % in inches
figure_height = 9; % in inches
screen_ppi = 72; 

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(fig, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);


%% HFC-134a and CFC-113/114/HCFC133a plot
subplot_height = 1/5.5;

% reported feedstock production
% prod_f_CFC113_1990_2022 = please contact UNEP Ozone Secretariat for reported feedstock production data
% prod_f_CFC114_1990_2022 = please contact UNEP Ozone Secretariat for reported feedstock production data

% total HFC-134a production by each pathway
pathway_1_HFC134a_nonA5 = pathway_fraction_nonA5(ResampleIndex_all_joint,:).*prod_HFC134a_nonA5;
pathway_2_HFC134a_nonA5 = (1-pathway_fraction_nonA5(ResampleIndex_all_joint,:)).*prod_HFC134a_nonA5;
pathway_1_HFC134a_A5 = pathway_fraction_A5(ResampleIndex_all_joint,:).*prod_HFC134a_A5;
pathway_2_HFC134a_A5 = (1-pathway_fraction_A5(ResampleIndex_all_joint,:)).*prod_HFC134a_A5;
pathway_1_HFC134a = pathway_1_HFC134a_A5 + pathway_1_HFC134a_nonA5;
pathway_2_HFC134a = pathway_2_HFC134a_A5 + pathway_2_HFC134a_nonA5;

vspace = 0.16;
hspace = 0.06;

fig = figure;
% top plot: HFC-134a production by each pathway
s1 = subplot_tight(4,8,3:6,[0,hspace]);
p1 = plot(1990:2022, prod_HFC134a_A5+prod_HFC134a_nonA5,'k','LineWidth',2);
LB = prctile(pathway_1_HFC134a,16);
MED = prctile(pathway_1_HFC134a,50);
UB  = prctile(pathway_1_HFC134a,84);
p2 = boundedline(1990:2022, MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0.3, 0.3, 0.3]);%[0.4660 0.6740 0.1880]);
LB = prctile(pathway_2_HFC134a,16);
MED = prctile(pathway_2_HFC134a,50);
UB  = prctile(pathway_2_HFC134a,84);
p3 = boundedline(1990:2022, MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0.7, 0.7, 0.7]);%[0.3010 0.7450 0.9330]);
title('A')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
ylabel({'HFC-134a','Production (Gg\cdoty^{-1})'})
% xlabel('Year')
xlim([2004,2019])
ylim([0,350])
leg = legend([p1,p3,p2],'Total','TCE pathway','PCE pathway', 'box' , 'on','edgecolor','w','NumColumns',1);
leg.ItemTokenSize(1) = 18;
legend boxoff 
set(gca,'FontSize',12)
legend('Location','eastoutside')
set(legend,'FontSize',12)
s1.Position(4) = subplot_height;

% second row: CFC-113, CFC-114, and HCFC-133a production
s4 = subplot_tight(3,8,10:11,[vspace,hspace]);
% p0 = plot(1990:2022,prod_f_CFC113_1990_2022,'Color','k','LineWidth',2);
LB = prctile(A5_113_FS_prod(ResampleIndex_CFC113,:),16);
MED = prctile(A5_113_FS_prod(ResampleIndex_CFC113,:),50);
UB = prctile(A5_113_FS_prod(ResampleIndex_CFC113,:),84);
p3 = boundedline(1990:2022,MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0.8500 0.3250 0.0980]);
hold on
LB = prctile(nonA5_113_FS_prod(ResampleIndex_CFC113,:),16);
MED = prctile(nonA5_113_FS_prod(ResampleIndex_CFC113,:),50);
UB = prctile(nonA5_113_FS_prod(ResampleIndex_CFC113,:),84);
p2 = boundedline(1990:2022,MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0 0 0.7]);
LB = prctile(total_113_FS_prod(ResampleIndex_CFC113,:),16);
MED = prctile(total_113_FS_prod(ResampleIndex_CFC113,:),50);
UB  = prctile(total_113_FS_prod(ResampleIndex_CFC113,:),84);
p1 = boundedline(1990:2022,MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0.3 0.3 0.3]);
xlim([2004,2019])
xticks([])
ylim([0,350])
yticks([0,100,200,300])
title('B    CFC-113')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
ylabel('Production (Gg\cdoty^{-1})')
% xlabel('Year')
set(gca,'FontSize',12)
s4.Position(4) = subplot_height;

% CFC-114 production
s3 = subplot_tight(3,8,12:13,[vspace,hspace]);
% p0 = plot(1990:2022,prod_f_CFC114_1990_2022,'Color','k','LineWidth',2);
LB = prctile(A5_114_FS_prod(ResampleIndex_all_joint,:),16);
MED = prctile(A5_114_FS_prod(ResampleIndex_all_joint,:),50);
UB = prctile(A5_114_FS_prod(ResampleIndex_all_joint,:),84);
p3 = boundedline(1990:2022,MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0.8500 0.3250 0.0980]);
hold on
LB = prctile(nonA5_114_FS_prod(ResampleIndex_all_joint,:),16);
MED = prctile(nonA5_114_FS_prod(ResampleIndex_all_joint,:),50);
UB = prctile(nonA5_114_FS_prod(ResampleIndex_all_joint,:),84);
p2 = boundedline(1990:2022,MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0 0 0.7]);
LB = prctile(total_114_FS_prod(ResampleIndex_all_joint,:),16);
MED = prctile(total_114_FS_prod(ResampleIndex_all_joint,:),50);
UB  = prctile(total_114_FS_prod(ResampleIndex_all_joint,:),84);
p1 = boundedline(1990:2022,MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0.3 0.3 0.3]);
xlim([2004,2019])
xticks([])
ylim([0,350])
% yticks([])
title("C    CFC-114")
ax = gca;
ax.TitleHorizontalAlignment = 'left';
% xlabel('Year')
set(gca,'FontSize',12)
s3.Position(4) = subplot_height;

% HCFC-133a production
total_133a = HCFC133a_HFC134a_A5 + HCFC133a_HFC134a_nonA5;
s4 = subplot_tight(3,8,14:15,[vspace,hspace]);
LB = prctile(HCFC133a_HFC134a_A5(ResampleIndex_all_joint,:),16);
MED = prctile(HCFC133a_HFC134a_A5(ResampleIndex_all_joint,:),50);
UB = prctile(HCFC133a_HFC134a_A5(ResampleIndex_all_joint,:),84);
p3 = boundedline(1990:2022,MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0.8500 0.3250 0.0980]);
hold on
LB = prctile(HCFC133a_HFC134a_nonA5(ResampleIndex_all_joint,:),16);
MED = prctile(HCFC133a_HFC134a_nonA5(ResampleIndex_all_joint,:),50);
UB = prctile(HCFC133a_HFC134a_nonA5(ResampleIndex_all_joint,:),84);
p2 = boundedline(1990:2022,MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0 0 0.7]);
LB = prctile(total_133a(ResampleIndex_all_joint,:),16);
MED = prctile(total_133a(ResampleIndex_all_joint,:),50);
UB  = prctile(total_133a(ResampleIndex_all_joint,:),84);
p1 = boundedline(1990:2022,MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0.3 0.3 0.3]);
xlim([2004,2019])
xticks([])
ylim([0,350])
% yticks([])
% title({'HCFC-133a','d) Production'})
title("D  HCFC-133a")
ax = gca;
ax.TitleHorizontalAlignment = 'left';
box on
% xlabel('Year')
set(gca,'FontSize',12)
leg = legend([p1,p2,p3],{'Global','Non-A 5','A 5'},'NumColumns',1);
% leg = legend([p1,p2,p3,p0],{'Global','Non-A 5','A 5',"Reported" + newline + "Feedstocks"},'NumColumns',1);
leg.ItemTokenSize(1) = 18;
legend boxoff 
legend('Location','eastoutside')
set(legend,'FontSize',12)
s4.Position(4) = subplot_height;

% bottom row: CFC-113, CFC-114, and HCFC-133a emissions
e_FS_CFC113 = e_A5_FS_CFC113 + e_nonA5_FS_CFC113;
s8 = subplot_tight(3,8,18:19,[vspace,hspace]);
p1 = plot(1982:2021,obs_derived_emiss_CFC113,'k','LineWidth',2,'LineStyle','--');
hold on
MED = prctile(e_A5_FS_CFC113(ResampleIndex_CFC113,1990-1954:2021-1954),50);
LB = prctile(e_A5_FS_CFC113(ResampleIndex_CFC113,1990-1954:2021-1954),16);
UB = prctile(e_A5_FS_CFC113(ResampleIndex_CFC113,1990-1954:2021-1954),84);
p5 = boundedline(1990:2021, MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.8500 0.3250 0.0980]);
MED = prctile(e_nonA5_FS_CFC113(ResampleIndex_CFC113,1990-1954:2021-1954),50);
LB = prctile(e_nonA5_FS_CFC113(ResampleIndex_CFC113,1990-1954:2021-1954),16);
UB = prctile(e_nonA5_FS_CFC113(ResampleIndex_CFC113,1990-1954:2021-1954),84);
p4 = boundedline(1990:2021, MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0. 0. 0.7]);
MED = prctile(e_FS_CFC113(ResampleIndex_CFC113,1990-1954:2021-1954),50);
LB = prctile(e_FS_CFC113(ResampleIndex_CFC113,1990-1954:2021-1954),16);
UB = prctile(e_FS_CFC113(ResampleIndex_CFC113,1990-1954:2021-1954),84);
p2 = boundedline(1990:2021, MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.3 0.3 0.3]);
xlim([2004,2019])
ylim([0,11])
yticks([0,5,10])
% xticks([])
% title('h) Emissions')
title("E")
ax = gca;
ax.TitleHorizontalAlignment = 'left';
ylabel('Emissions (Gg\cdoty^{-1})')
% xlabel('Year')
set(gca,'FontSize',12)
s8.Position(4) = subplot_height;

e_FS_CFC114 = e_A5_FS_CFC114 + e_nonA5_FS_CFC114;
s9 = subplot_tight(3,8,20:21,[vspace,hspace]);
p1 = plot(1978:2021,obs_derived_emiss_CFC114(1:44),'k','LineWidth',2,'LineStyle','--');
hold on
MED = prctile(e_A5_FS_CFC114(ResampleIndex_all_joint,1990-1934:2021-1934),50);
LB = prctile(e_A5_FS_CFC114(ResampleIndex_all_joint,1990-1934:2021-1934),16);
UB = prctile(e_A5_FS_CFC114(ResampleIndex_all_joint,1990-1934:2021-1934),84);
p5 = boundedline(1990:2021, MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.8500 0.3250 0.0980]);
MED = prctile(e_nonA5_FS_CFC114(ResampleIndex_all_joint,1990-1934:2021-1934),50);
LB = prctile(e_nonA5_FS_CFC114(ResampleIndex_all_joint,1990-1934:2021-1934),16);
UB = prctile(e_nonA5_FS_CFC114(ResampleIndex_all_joint,1990-1934:2021-1934),84);
p4 = boundedline(1990:2021, MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0. 0. 0.7]);
MED = prctile(e_FS_CFC114(ResampleIndex_all_joint,1990-1934:2021-1934),50);
LB = prctile(e_FS_CFC114(ResampleIndex_all_joint,1990-1934:2021-1934),16);
UB = prctile(e_FS_CFC114(ResampleIndex_all_joint,1990-1934:2021-1934),84);
p2 = boundedline(1990:2021, MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.3 0.3 0.3]);
xlim([2004,2019])
ylim([0,3.5])
% xticks([])
% title('i) Emissions')
title("F")
ax = gca;
ax.TitleHorizontalAlignment = 'left';
xlabel('Year')
set(gca,'FontSize',12)
s9.Position(4) = subplot_height;

e_HFC134a_HCFC133a = e_HFC134a_HCFC133a_A5 + e_HFC134a_HCFC133a_nonA5;
s10 = subplot_tight(3,8,22:23,[vspace,hspace]);
p1 = plot(1990:2021,obs_derived_emiss_HCFC133a,'k','LineWidth',2,'LineStyle','--');
hold on
MED = prctile(e_HFC134a_HCFC133a_A5(ResampleIndex_all_joint,1:33),50);
LB = prctile(e_HFC134a_HCFC133a_A5(ResampleIndex_all_joint,1:33),16);
UB = prctile(e_HFC134a_HCFC133a_A5(ResampleIndex_all_joint,1:33),84);
p5 = boundedline(1990:2022, MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.8500 0.3250 0.0980]);
MED = prctile(e_HFC134a_HCFC133a_nonA5(ResampleIndex_all_joint,1:33),50);
LB = prctile(e_HFC134a_HCFC133a_nonA5(ResampleIndex_all_joint,1:33),16);
UB = prctile(e_HFC134a_HCFC133a_nonA5(ResampleIndex_all_joint,1:33),84);
p4 = boundedline(1990:2022, MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0. 0. 0.7]);
MED = prctile(e_HFC134a_HCFC133a(ResampleIndex_all_joint,1:33),50);
LB = prctile(e_HFC134a_HCFC133a(ResampleIndex_all_joint,1:33),16);
UB = prctile(e_HFC134a_HCFC133a(ResampleIndex_all_joint,1:33),84);
p2 = boundedline(1990:2022, MED',[MED'- LB',UB' - MED'],'alpha','cmap',[0.3 0.3 0.3]);
xlim([2004,2019])
ylim([0,3.8])
% xticks([])
% title('j) Emissions')
title("G")
ax = gca;
ax.TitleHorizontalAlignment = 'left';
% xlabel('Year')
leg = legend([p1,p2,p4,p5],"Obs." + newline + "derived",'Global','Non-A 5','A 5');
leg.ItemTokenSize(1) = 18;
legend boxoff 
legend('Location','eastoutside')
set(legend,'FontSize',12)
set(gca,'FontSize',12)
s10.Position(4) = subplot_height;

figure_width = 8; % in inches
figure_height = 6; % in inches
screen_ppi = 72; 

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(fig, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

%% HFC-125 and CFC-113/114/115 plot

% Plot emission rates and emissions
vspace = 0.05;
hspace = 0.05;

fig = figure;
subplot_height = 1/3;

% CFC-113 emissions
subplot(2,2,1)
p1 = plot(1982:2021,obs_derived_emiss_CFC113,'k','LineWidth',2,'LineStyle','--');
LB = prctile(e_HFC125_CFC113(ResampleIndex_CFC113,:),16);
MED = prctile(e_HFC125_CFC113(ResampleIndex_CFC113,:),50);
UB  = prctile(e_HFC125_CFC113(ResampleIndex_CFC113,:),84);
p2 = boundedline(1955:2020,MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0. 0. 0.]);
hold on
LB = prctile(e_HFC125_CFC113_nonA5(ResampleIndex_CFC113,:),16);
MED = prctile(e_HFC125_CFC113_nonA5(ResampleIndex_CFC113,:),50);
UB = prctile(e_HFC125_CFC113_nonA5(ResampleIndex_CFC113,:),84);
p4 = boundedline(1990:2022,MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0 0 0.7]);
LB = prctile(e_HFC125_CFC113_A5(ResampleIndex_CFC113,:),16);
MED = prctile(e_HFC125_CFC113_A5(ResampleIndex_CFC113,:),50);
UB = prctile(e_HFC125_CFC113_A5(ResampleIndex_CFC113,:),84);
p5 = boundedline(1990:2022,MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0.8500 0.3250 0.0980]);
xlim([2004,2019])
xticks([])
yticks([0,4,8])
ylim([0,9])
ylabel('Emissions (Gg\cdoty^{-1})')
title('A               CFC-113')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
set(gca,'FontSize',12)

% CFC-114 emissions
subplot(2,2,2)
p1 = plot(1978:2021,obs_derived_emiss_CFC114(1:44),'k','LineWidth',2,'LineStyle','--');
LB = prctile(e_HFC125_CFC114(ResampleIndex_all_joint,:),16);
MED = prctile(e_HFC125_CFC114(ResampleIndex_all_joint,:),50);
UB  = prctile(e_HFC125_CFC114(ResampleIndex_all_joint,:),84);
p2 = boundedline(1935:2020,MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0.3 0.3 0.3]);
hold on
LB = prctile(e_HFC125_CFC114_nonA5(ResampleIndex_all_joint,:),16);
MED = prctile(e_HFC125_CFC114_nonA5(ResampleIndex_all_joint,:),50);
UB = prctile(e_HFC125_CFC114_nonA5(ResampleIndex_all_joint,:),84);
p4 = boundedline(1990:2022,MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0 0 0.7]);
LB = prctile(e_HFC125_CFC114_A5(ResampleIndex_all_joint,:),16);
MED = prctile(e_HFC125_CFC114_A5(ResampleIndex_all_joint,:),50);
UB = prctile(e_HFC125_CFC114_A5(ResampleIndex_all_joint,:),84);
p5 = boundedline(1990:2022,MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0.8500 0.3250 0.0980]);
xlim([2004,2019])
xticks([])
ylim([0,3.6])
title('B               CFC-114')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
set(gca,'FontSize',12)

% CFC-115 emissions
subplot(2,2,3)
p1 = plot(1978:2021,obs_derived_emiss_CFC115(1:44),'k','LineWidth',2,'LineStyle','--');
LB = prctile(e_HFC125_CFC115(ResampleIndex_CFC115,:),16);
MED = prctile(e_HFC125_CFC115(ResampleIndex_CFC115,:),50);
UB  = prctile(e_HFC125_CFC115(ResampleIndex_CFC115,:),84);
p2 = boundedline(1935:2020,MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0.3 0.3 0.3]);
hold on
LB = prctile(e_HFC125_CFC115_nonA5(ResampleIndex_CFC115,:),16);
MED = prctile(e_HFC125_CFC115_nonA5(ResampleIndex_CFC115,:),50);
UB = prctile(e_HFC125_CFC115_nonA5(ResampleIndex_CFC115,:),84);
p4 = boundedline(1990:2022,MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0 0 0.7]);
LB = prctile(e_HFC125_CFC115_A5(ResampleIndex_CFC115,:),16);
MED = prctile(e_HFC125_CFC115_A5(ResampleIndex_CFC115,:),50);
UB = prctile(e_HFC125_CFC115_A5(ResampleIndex_CFC115,:),84);
p5 = boundedline(1990:2022,MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0.8500 0.3250 0.0980]);
xlim([2004,2019])
ylim([0,2.4])
yticks([0,1,2])
ylabel('Emissions (Gg\cdoty^{-1})')
xlabel('Year')
title('C               CFC-115')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
set(gca,'FontSize',12)
box on

% HCFC-133a emissions
e_HFC125_HCFC133a = e_HFC125_HCFC133a_A5 + e_HFC125_HCFC133a_nonA5;
subplot(2,2,4)
p1 = plot(1990:2021,obs_derived_emiss_HCFC133a,'k','LineWidth',2,'LineStyle','--');
hold on
LB = prctile(e_HFC125_HCFC133a(ResampleIndex_all_joint,:),16);
MED = prctile(e_HFC125_HCFC133a(ResampleIndex_all_joint,:),50);
UB  = prctile(e_HFC125_HCFC133a(ResampleIndex_all_joint,:),84);
p2 = boundedline(1990:2022,MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0.3 0.3 0.3]);
LB = prctile(e_HFC125_HCFC133a_nonA5(ResampleIndex_all_joint,:),16);
MED = prctile(e_HFC125_HCFC133a_nonA5(ResampleIndex_all_joint,:),50);
UB = prctile(e_HFC125_HCFC133a_nonA5(ResampleIndex_all_joint,:),84);
p4 = boundedline(1990:2022,MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0 0 0.7]);
LB = prctile(e_HFC125_HCFC133a_A5(ResampleIndex_all_joint,:),16);
MED = prctile(e_HFC125_HCFC133a_A5(ResampleIndex_all_joint,:),50);
UB = prctile(e_HFC125_HCFC133a_A5(ResampleIndex_all_joint,:),84);
p5 = boundedline(1990:2022,MED',[MED'-LB',UB'-MED'],'alpha','cmap',[0.8500 0.3250 0.0980]);
xlim([2004,2019])
ylim([0,3.8])
% title('i) Emissions')
title('D             HCFC-133a')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
xlabel('Year')
% ylabel('Gg/yr')
leg = legend([p1,p2,p4,p5],{'Obs. derived','Global','Non-A 5','A 5'},'NumColumns',1);
leg.ItemTokenSize(1) = 18;
legend boxoff 
legend('Location','eastoutside')
set(legend,'FontSize',12)
set(gca,'FontSize',12)
box on

figure_width = 8; % in inches
figure_height = 6; % in inches
screen_ppi = 72; 

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(fig, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

%% fraction of CFC-113 and CFC-114 produced from HFC-125
% HFC125_fraction_CFC113 = sum(mean(e_HFC125_CFC113_constrained(ResampleIndex_CFC113_constrained,end-18:end)))...
% /sum(mean(emissions_CFC113_constrained(ResampleIndex_CFC113_constrained,end-18:end)))
% 
% HFC125_fraction_CFC114 = sum(mean(e_HFC125_CFC114(ResampleIndex_all_joint,end-18:end)))...
% /sum(mean(emissions_CFC114_total(ResampleIndex_all_joint,end-18:end)))
% 

%% ODP and CO2eq as timeseries
% ODP and GWP from https://ozone.unep.org/sites/default/files/2023-02/Scientific-Assessment-of-Ozone-Depletion-2022.pdf

% mass of CFC-11 for conversion
CFC11_mass = 137.37;

% ODP from CFC113 and 113a
ODP_CFC113 = 0.78;
ODP_CFC113_HFC134a = ODP_CFC113*CFC11_mass/CFC113_mass*e_FS_CFC113(ResampleIndex_CFC113,:); 
ODP_CFC113_HFC125 = ODP_CFC113*CFC11_mass/CFC113_mass*e_HFC125_CFC113(ResampleIndex_CFC113,:); 

% ODP from CFC114
ODP_CFC114 = 0.61;
ODP_CFC114_HFC134a = ODP_CFC114*CFC11_mass/CFC114_mass*e_FS_CFC114(ResampleIndex_all_joint,:); 
ODP_CFC114_HFC125 = ODP_CFC114*CFC11_mass/CFC114_mass*e_HFC125_CFC114(ResampleIndex_all_joint,:); 

% ODP from CFC115
ODP_CFC115 = 0.45;
ODP_CFC115_HFC125 = ODP_CFC115*CFC11_mass/CFC115_mass*e_HFC125_CFC115(ResampleIndex_CFC115,:); 

% ODP from HCFC133a
ODP_HCFC133a = 0.019;
ODP_HCFC133a_HFC134a = ODP_HCFC133a*CFC11_mass/HCFC133a_mass*e_HFC134a_HCFC133a(ResampleIndex_all_joint,:); 
ODP_HCFC133a_HFC125 = ODP_HCFC133a*CFC11_mass/HCFC133a_mass*e_HFC125_HCFC133a(ResampleIndex_all_joint,:); 

% ODP for all species combined
ODP_HFC134a = ODP_CFC113_HFC134a(:,36:end)+ODP_CFC114_HFC134a(:,56:end)+ODP_HCFC133a_HFC134a;
ODP_HFC125 = ODP_CFC113_HFC125(:,36:end)+ODP_CFC114_HFC125(:,56:end)+ODP_CFC115_HFC125(1:Nresamps2,56:end)+ODP_HCFC133a_HFC125(:,1:end-2);

% GWP
%CFC-113
GWP_CFC113 = 5490;
GWP_CFC113_HFC134a = GWP_CFC113*e_FS_CFC113(ResampleIndex_CFC113,:); 
GWP_CFC113_HFC125 = GWP_CFC113*e_HFC125_CFC113(ResampleIndex_CFC113,:); 

%CFC-114
GWP_CFC114 = 8634;
GWP_CFC114_HFC134a = GWP_CFC114*e_FS_CFC114(ResampleIndex_all_joint,:); 
GWP_CFC114_HFC125 = GWP_CFC114*e_HFC125_CFC114(ResampleIndex_all_joint,:); 

%CFC-115
GWP_CFC115 = 9630;
GWP_CFC115_HFC125 = GWP_CFC115*e_HFC125_CFC115(ResampleIndex_CFC115,:); 

%HCFC-133a
GWP_HCFC133a = 378;
GWP_HCFC133a_HFC134a = GWP_HCFC133a*e_HFC134a_HCFC133a(ResampleIndex_all_joint,:); 
GWP_HCFC133a_HFC125 = GWP_HCFC133a*e_HFC125_HCFC133a(ResampleIndex_all_joint,:); 

% GWP for all species combined
GWP_HFC134a = GWP_CFC113_HFC134a(:,36:end)+GWP_CFC114_HFC134a(:,56:end)+GWP_HCFC133a_HFC134a;
GWP_HFC125 = GWP_CFC113_HFC125(:,36:end)+GWP_CFC114_HFC125(:,56:end)+GWP_CFC115_HFC125(1:Nresamps2,56:end)+GWP_HCFC133a_HFC125(:,1:end-2);

% GWP of HFC emissions
% GWP of HFC-134a: 1330
% GWP of HFC-125: 3820
e_HFC134a = 1/1000*[1415.8,3599.3,5578.9,7700.1,13770,21655.3,33888.6,40358.2,59579.2,67798.6,74935,81770.9,93533.7,103006.6, ...
    110787.8,114222.7,118512.8,129036.5,136277.8,143907.5,155912.4,157934.4,164357.9,171714.9,183725.3,196012.6,211348.9, ...
    219816,219957.8,223178.1];
e_HFC_125 = 1/1000*[2.4,8.2,615.6,1037,1069.9,3111,4213.9,7576.5,6231.5,5881.8,8222.6,9201.6,10844.9,12917.9,14031.7,15333.6,18438.2, ...
    21449.2,24836.7,27617.5,33607.6,38192.2,43229.4,48158.7,54146.2,58970.4,65888.1,74919.3,79151.5,82654];
GWP_HFC134a_1990_2019 = 1300*e_HFC134a;
GWP_HFC125_1990_2019 = 3170*e_HFC_125;

%% Plot of ODP and GWP
fig = figure;
subplot(2,2,1)
p1 = area(1990:2022,prctile(ODP_HFC134a,50),'FaceColor',[0.4940 0.1840 0.5560]);
hold on
p2 = area(1990:2022,prctile(ODP_CFC113_HFC134a(:,36:end)+ODP_CFC114_HFC134a(:,56:end),50),'FaceColor',[0.8500 0.3250 0.0980]);
p3 = area(1990:2022,prctile(ODP_CFC113_HFC134a(:,36:end),50),'FaceColor',[0 0.4470 0.7410]);
leg = legend([p3,p2,p1],'CFC-113','CFC-114','HCFC-133a','Location','northwest');
xlim([2004,2019])
ylim([0,4.2])
title('A             HFC-134a')
ax = gca;
ax.TitleHorizontalAlignment = 'left';
ylabel({'CFC-11 Equivalent', 'Emissions (ODP-Gg\cdotyr^{-1})'})
leg.ItemTokenSize(1) = 18;
legend boxoff 
set(gca,'FontSize',12)
set(legend,'FontSize',12)
set(gca,'FontSize',12)

subplot(2,2,2)
p1 = area(1990:2020,prctile(ODP_HFC125,50),'FaceColor',[0.4940 0.1840 0.5560]);
hold on
p2 = area(1990:2020,prctile(ODP_CFC113_HFC125(:,36:end)+ODP_CFC114_HFC125(:,56:end)+ODP_CFC115_HFC125(1:Nresamps2,56:end),50),'FaceColor',[0.4660 0.6740 0.1880]);
p3 = area(1990:2020,prctile(ODP_CFC113_HFC125(:,36:end)+ODP_CFC114_HFC125(:,56:end),50),'FaceColor',[0.8500 0.3250 0.0980]);
p4 = area(1990:2020,prctile(ODP_CFC113_HFC125(:,36:end),50),'FaceColor',[0 0.4470 0.7410]);
leg = legend([p4,p3,p2,p1],'CFC-113','CFC-114','CFC-115','HCFC-133a','Location','northwest');
xlim([2004,2019])
ylim([0,1.2])
leg.ItemTokenSize(1) = 18;
legend boxoff 
set(gca,'FontSize',12)
set(legend,'FontSize',12)
set(gca,'FontSize',12)
title('B                  HFC-125')
ax = gca;
ax.TitleHorizontalAlignment = 'left';

subplot(2,2,3)
p1 = area(1990:2019,GWP_HFC134a_1990_2019/1000+prctile(GWP_HFC134a(:,1:end-3),50)/1000,'FaceColor',[0.4940 0.1840 0.5560]);
hold on
p2 = area(1990:2019,GWP_HFC134a_1990_2019/1000+prctile(GWP_CFC113_HFC134a(:,36:end-3)+GWP_CFC114_HFC134a(:,56:end-3),50)/1000,'FaceColor',[0.8500 0.3250 0.0980]);
p3 = area(1990:2019,GWP_HFC134a_1990_2019/1000+prctile(GWP_CFC113_HFC134a(:,36:end-3),50)/1000,'FaceColor',[0 0.4470 0.7410]);
p4 = area(1990:2019,GWP_HFC134a_1990_2019/1000,'FaceColor',[0.7,0.7,0.7]);
xlim([2004,2019])
ylim([0, 360])
leg = legend([p4,p3,p2,p1],'HFC-134a','CFC-113','CFC-114','HCFC-133a','Location','northwest');
leg.ItemTokenSize(1) = 18;
legend boxon 
set(gca,'FontSize',12)
set(legend,'FontSize',12)
set(gca,'FontSize',12)
xlabel('Year')
ylabel({'CO_2 Equivalent', 'Emissions (TgCO_2e\cdotyr^{-1})'})
title('C')
ax = gca;
ax.TitleHorizontalAlignment = 'left';

subplot(2,2,4)
p1 = area(1990:2019,GWP_HFC125_1990_2019/1000+prctile(GWP_HFC125(:,1:end-1),50)/1000,'FaceColor',[0.4940 0.1840 0.5560]);
hold on
p2 = area(1990:2019,GWP_HFC125_1990_2019/1000+prctile(GWP_CFC113_HFC125(:,36:end-1)+GWP_CFC114_HFC125(:,56:end-1)+GWP_CFC115_HFC125(1:Nresamps2,56:end-1),50)/1000,'FaceColor',[0.4660 0.6740 0.1880]);
p3 = area(1990:2019,GWP_HFC125_1990_2019/1000+prctile(GWP_CFC113_HFC125(:,36:end-1)+GWP_CFC114_HFC125(:,56:end-1),50)/1000,'FaceColor',[0.8500 0.3250 0.0980]);
p4 = area(1990:2019,GWP_HFC125_1990_2019/1000+prctile(GWP_CFC113_HFC125(:,36:end-1),50)/1000,'FaceColor',[0 0.4470 0.7410]);
p5 = area(1990:2019,GWP_HFC125_1990_2019/1000,'FaceColor',[0.7,0.7,0.7]);
xlim([2004,2019])
ylim([0,320])
xlabel('Year')
leg = legend([p5,p4,p3,p2,p1],'HFC-125','CFC-113','CFC-114','CFC-115','HCFC-133a');
leg.ItemTokenSize(1) = 18;
legend boxoff 
set(gca,'FontSize',12)
set(legend,'FontSize',12)
set(gca,'FontSize',12)
title('D')
ax = gca;
ax.TitleHorizontalAlignment = 'left';

figure_width = 8; % in inches
figure_height = 6; % in inches
screen_ppi = 72; 

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(fig, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);


%% future GWP and ODP (projection from Velders et al 2022)

HFC_134a_2020_2050_LB = 1/1000*[323799.1,330883.1,342132.7,356170.4,341415,351976.4, ...
    359269,364359.6,368263.5,325682.8,293627,296898.1,294826.7,296606.1,290864.5, ...
    239369.3,235447.3,232496,233669.9,235203.3,178430.8,178226.9,174090.8,173899.3, ...
    173709.9,91934.5,91850.3,70218.9,70140.7,70061,69979.5];
HFC_134a_2020_2050_UB = 1/1000*[317332.7,318537.5,324383.8,332373,310723.2,325325.6, ...
    331494.8,335340.1,340949.5,307288.8,293169.2,295490.2,294529.3,296258.1,290333.4, ...
    239062.6,235161.3,231890.5,232887.3,233895.3,179780.7,180405.2,176660.3,177248.8, ...
    177839.6,96124.7,96334,72326.4,72473.2,72609.3,72745.3];
HFC_134a_2020_2050 = 0.5*(HFC_134a_2020_2050_LB + HFC_134a_2020_2050_UB);

HFC_125_2020_2050_LB = 1/1000*[205100.4,209308.2,217918.3,227000.9,215292.4,251908, ...
    253991.2,254316.3,256642.3,229213.5,242829.4,243318.6,239086.6,238223.2,230574.6, ...
    193467.8,188559.6,182368.7,181930.9,181494.8,145564.6,145286.6,139292.5,139024.6, ...
    138754.7,85259.2,85147.9,53687,53598.1,53511.3,53423.6];
HFC_125_2020_2050_UB = 1/1000*[216789.2,228117.6,242945.7,254909.3,239590.3,287315.1, ...
    290083.4,291140.9,288542.8,252528.5,254714.8,252682.4,244210.6,243021.8,235025.8, ...
    198999.1,193766.2,186056.1,185283.3,184830.8,150653.7,150603.5,143627.5,143562.9, ...
    143491.6,91752.5,91627.7,53465.6,53348.8,53226.6,53098.8];
HFC_125_2020_2050 = 0.5*(HFC_125_2020_2050_LB + HFC_125_2020_2050_UB);

% vary emission rates to vary totals 
% medians: 0.04, 0.015, 0.006
% 16%: 0.026, 0.009, 0.004
% 84%: 0.054, 0.021, 0.01
CFC113_HFC134a_future = 0.04*(1-0.65)*HFC_134a_2020_2050;
CFC114_HFC134a_future = 0.015*(1-0.65)*HFC_134a_2020_2050;
HCFC133a_HFC134a_future = 0.006*(0.65)*HFC_134a_2020_2050;

CFC113_HFC134a_future_ODP = ODP_CFC113*CFC11_mass/CFC113_mass*CFC113_HFC134a_future;
CFC114_HFC134a_future_ODP = ODP_CFC114*CFC11_mass/CFC114_mass*CFC114_HFC134a_future;
HCFC133a_HFC134a_future_ODP = ODP_HCFC133a*CFC11_mass/HCFC133a_mass*HCFC133a_HFC134a_future;
CFC113_HFC134a_future_GWP = GWP_CFC113*CFC113_HFC134a_future;
CFC114_HFC134a_future_GWP = GWP_CFC114*CFC114_HFC134a_future;
HCFC133a_HFC134a_future_GWP = GWP_HCFC133a*HCFC133a_HFC134a_future;

HFC134a_future_ODP = CFC113_HFC134a_future_ODP+CFC114_HFC134a_future_ODP+HCFC133a_HFC134a_future_ODP;
HFC134a_future_GWP = CFC113_HFC134a_future_GWP+CFC114_HFC134a_future_GWP+HCFC133a_HFC134a_future_GWP;

% vary emission rates to vary totals
% medians: 0.002, 0.002, 0.007, 0.007
% 16%: 0.001, 0.001, 0.005, 0.004
% 84%: 0.003, 0.003, 0.010, 0.009
CFC113_HFC125_future = 0.002*HFC_125_2020_2050;
CFC114_HFC125_future = 0.002*HFC_125_2020_2050;
CFC115_HFC125_future = 0.007*HFC_125_2020_2050;
HCFC133a_HFC125_future = 0.007*HFC_125_2020_2050;

CFC113_HFC125_future_ODP = ODP_CFC113*CFC11_mass/CFC113_mass*CFC113_HFC125_future;
CFC114_HFC125_future_ODP = ODP_CFC114*CFC11_mass/CFC114_mass*CFC114_HFC125_future;
CFC115_HFC125_future_ODP = ODP_CFC115*CFC11_mass/CFC115_mass*CFC115_HFC125_future;
HCFC133a_HFC125_future_ODP = ODP_HCFC133a*CFC11_mass/HCFC133a_mass*HCFC133a_HFC125_future;
CFC113_HFC125_future_GWP = GWP_CFC113*CFC113_HFC125_future;
CFC114_HFC125_future_GWP = GWP_CFC114*CFC114_HFC125_future;
CFC115_HFC125_future_GWP = GWP_CFC114*CFC115_HFC125_future;
HCFC133a_HFC125_future_GWP = GWP_HCFC133a*HCFC133a_HFC125_future;

HFC125_future_ODP = CFC113_HFC125_future_ODP+CFC114_HFC125_future_ODP+CFC115_HFC125_future_ODP+HCFC133a_HFC125_future_ODP;
HFC125_future_GWP = CFC113_HFC125_future_GWP+CFC114_HFC125_future_GWP+CFC115_HFC125_future_GWP+HCFC133a_HFC125_future_GWP;

sum(HFC125_future_ODP)
sum(HFC125_future_GWP/1000)
sum(HFC134a_future_ODP)
sum(HFC134a_future_GWP/1000)

%% ODP and CO2eq per Gg produced
% ODP
ODP_CFC113_HFC134a = ODP_CFC113*e_FS_CFC113(ResampleIndex_CFC113,end-32:end)./(prod_HFC134a_A5+prod_HFC134a_nonA5);
% CFC-114
ODP_CFC114_HFC134a = ODP_CFC114*e_FS_CFC114(ResampleIndex_all_joint,end-32:end)./(prod_HFC134a_A5+prod_HFC134a_nonA5);
% HCFC-133a
ODP_HCFC133a_HFC134a = ODP_HCFC133a*e_HFC134a_HCFC133a(ResampleIndex_all_joint,:)./(prod_HFC134a_A5+prod_HFC134a_nonA5);
% total
ODP_HFC134a = ODP_CFC113_HFC134a+ODP_CFC114_HFC134a+ODP_HCFC133a_HFC134a;

% global warming potential per unit of HFC-125 produced
% CFC-113
ODP_CFC113_HFC125 = ODP_CFC113*e_HFC125_CFC113(ResampleIndex_CFC113,end-32:end)./(prod_HFC125_A5+prod_HFC125_nonA5);
% CFC-114
ODP_CFC114_HFC125 = ODP_CFC114*e_HFC125_CFC114(ResampleIndex_all_joint,end-32:end)./(prod_HFC125_A5+prod_HFC125_nonA5);
% CFC-115
ODP_CFC115_HFC125 = ODP_CFC115*e_HFC125_CFC115(ResampleIndex_CFC115,end-32:end)./(prod_HFC125_A5+prod_HFC125_nonA5);
% HCFC-133a
ODP_HCFC133a_HFC125 = ODP_HCFC133a*e_HFC125_HCFC133a(ResampleIndex_all_joint,:)./(prod_HFC125_A5+prod_HFC125_nonA5);
% total
ODP_HFC125 = ODP_CFC113_HFC125+ODP_CFC114_HFC125+ODP_CFC115_HFC125(1:Nresamps2,:)+ODP_HCFC133a_HFC125;


% global warming potential per unit of HFC-134a produced
% CFC-113
GWP_CFC113_HFC134a = GWP_CFC113*e_FS_CFC113(ResampleIndex_CFC113,end-32:end)./(prod_HFC134a_A5+prod_HFC134a_nonA5);
% CFC-114
GWP_CFC114_HFC134a = GWP_CFC114*e_FS_CFC114(ResampleIndex_all_joint,end-32:end)./(prod_HFC134a_A5+prod_HFC134a_nonA5);
% HCFC-133a
GWP_HCFC133a_HFC134a = GWP_HCFC133a*e_HFC134a_HCFC133a(ResampleIndex_all_joint,:)./(prod_HFC134a_A5+prod_HFC134a_nonA5);
% total
GWP_HFC134a = GWP_CFC113_HFC134a+GWP_CFC114_HFC134a+GWP_HCFC133a_HFC134a;

% global warming potential per unit of HFC-125 produced
% CFC-113
GWP_CFC113_HFC125 = GWP_CFC113*e_HFC125_CFC113(ResampleIndex_CFC113,end-32:end)./(prod_HFC125_A5+prod_HFC125_nonA5);
% CFC-114
GWP_CFC114_HFC125 = GWP_CFC114*e_HFC125_CFC114(ResampleIndex_all_joint,end-32:end)./(prod_HFC125_A5+prod_HFC125_nonA5);
% CFC-115
GWP_CFC115_HFC125 = GWP_CFC115*e_HFC125_CFC115(ResampleIndex_CFC115,end-32:end)./(prod_HFC125_A5+prod_HFC125_nonA5);
% HCFC-133a
GWP_HCFC133a_HFC125 = GWP_HCFC133a*e_HFC125_HCFC133a(ResampleIndex_all_joint,:)./(prod_HFC125_A5+prod_HFC125_nonA5);
% total
GWP_HFC125 = GWP_CFC113_HFC125+GWP_CFC114_HFC125+GWP_CFC115_HFC125(1:Nresamps2,:)+GWP_HCFC133a_HFC125;


%% emission rates relative to HFC production
e_HFC134a_CFC113_per_HFC134a = e_FS_CFC113(ResampleIndex_CFC113,end-5:end)./pathway_1_HFC134a(:,end-5:end);
e_HFC134a_CFC114_per_HFC134a = e_FS_CFC114(ResampleIndex_all_joint,end-5:end)./pathway_1_HFC134a(:,end-5:end);
e_HFC134a_HCFC133a_per_HFC134a = e_HFC134a_HCFC133a(ResampleIndex_all_joint,end-5:end)./pathway_2_HFC134a(:,end-5:end);

mean(prctile(e_HFC134a_CFC113_per_HFC134a,16)),mean(prctile(e_HFC134a_CFC113_per_HFC134a,50)),mean(prctile(e_HFC134a_CFC113_per_HFC134a,84))
mean(prctile(e_HFC134a_CFC114_per_HFC134a,16)),mean(prctile(e_HFC134a_CFC114_per_HFC134a,50)),mean(prctile(e_HFC134a_CFC114_per_HFC134a,84))
mean(prctile(e_HFC134a_HCFC133a_per_HFC134a,16)),mean(prctile(e_HFC134a_HCFC133a_per_HFC134a,50)),mean(prctile(e_HFC134a_HCFC133a_per_HFC134a,84))

%% emission rates percentiles
% global emission rates
emission_rate_CFC113_FS_global = e_FS_CFC113(ResampleIndex_CFC113,1990-1954:2022-1954)./total_113_FS_prod(ResampleIndex_CFC113,:);
emission_rate_CFC114_FS_global = e_FS_CFC114(ResampleIndex_all_joint,1990-1934:2022-1934)./total_114_FS_prod(ResampleIndex_all_joint,:);
emission_rate_HCFC133a_FS_global = e_HFC134a_HCFC133a(ResampleIndex_all_joint,:)./total_133a(ResampleIndex_all_joint,:);

% Global numbers reported in Table 1
mean(prctile(emission_rate_CFC113_FS_global(:,end-5:end),16)),mean(prctile(emission_rate_CFC114_FS_global(:,end-5:end),16)),mean(prctile(emission_rate_HCFC133a_FS_global(:,end-5:end),16))
mean(prctile(emission_rate_CFC113_FS_global(:,end-5:end),50)),mean(prctile(emission_rate_CFC114_FS_global(:,end-5:end),50)),mean(prctile(emission_rate_HCFC133a_FS_global(:,end-5:end),50))
mean(prctile(emission_rate_CFC113_FS_global(:,end-5:end),84)),mean(prctile(emission_rate_CFC114_FS_global(:,end-5:end),84)),mean(prctile(emission_rate_HCFC133a_FS_global(:,end-5:end),84))

% non-A 5 numbers in Table 1
prctile(emission_rate_CFC113_FS_nonA5(ResampleIndex_CFC113),16)
prctile(emission_rate_CFC113_FS_nonA5(ResampleIndex_CFC113),50)
prctile(emission_rate_CFC113_FS_nonA5(ResampleIndex_CFC113),84)
prctile(emission_rate_CFC114_FS_nonA5(ResampleIndex_all_joint),16)
prctile(emission_rate_CFC114_FS_nonA5(ResampleIndex_all_joint),50)
prctile(emission_rate_CFC114_FS_nonA5(ResampleIndex_all_joint),84)
prctile(emission_rate_HFC134a_HCFC133a_nonA5(ResampleIndex_all_joint),16)
prctile(emission_rate_HFC134a_HCFC133a_nonA5(ResampleIndex_all_joint),50)
prctile(emission_rate_HFC134a_HCFC133a_nonA5(ResampleIndex_all_joint),84)

% emission rates relative to HFC134a production
mean(prctile(e_FS_CFC113(ResampleIndex_CFC113,end-8:end-3)./pathway_1_HFC134a(:,end-8:end-3),50))
mean(prctile(e_FS_CFC114(ResampleIndex_all_joint,end-8:end-3)./pathway_1_HFC134a(:,end-8:end-3),50))
mean(prctile(e_HFC134a_HCFC133a(ResampleIndex_all_joint,end-8:end-3)./pathway_2_HFC134a(:,end-8:end-3),50))

% emission rates relative to HFC125 production, global
mean(prctile(e_HFC125_CFC113(ResampleIndex_CFC113,end-8:end-3)./(prod_HFC125_A5(:,end-8:end-3)+prod_HFC125_nonA5(:,end-8:end-3)),50))
mean(prctile(e_HFC125_CFC114(ResampleIndex_all_joint,end-8:end-3)./(prod_HFC125_A5(:,end-8:end-3)+prod_HFC125_nonA5(:,end-8:end-3)),50))
mean(prctile(e_HFC125_CFC115(ResampleIndex_CFC115,end-8:end-3)./(prod_HFC125_A5(:,end-8:end-3)+prod_HFC125_nonA5(:,end-8:end-3)),50))
mean(prctile(e_HFC125_HCFC133a(ResampleIndex_all_joint,end-8:end-3)./(prod_HFC125_A5(:,end-8:end-3)+prod_HFC125_nonA5(:,end-8:end-3)),50))

% emission rates relative to HFC125 production, non-A5
prctile(byproduct_rate_125_113_nonA5(ResampleIndex_CFC113,end),50)
prctile(byproduct_rate_125_114_nonA5(ResampleIndex_all_joint,end),50)
prctile(byproduct_rate_125_115_nonA5(ResampleIndex_CFC115,end),50)
prctile(emission_rate_HFC125_HCFC133a_nonA5(ResampleIndex_all_joint,end),50)

%% fraction of emissions by source
CFC_113_HFC134a_fraction = e_FS_CFC113(ResampleIndex_CFC113,1:end-3)./(emissions_CFC113(ResampleIndex_CFC113,1:end));
% 2015-2019
mean(prctile(CFC_113_HFC134a_fraction(:,end-4:end),16)),mean(prctile(CFC_113_HFC134a_fraction(:,end-4:end),50)),mean(prctile(CFC_113_HFC134a_fraction(:,end-4:end),84))
% 2004-2019
mean(prctile(CFC_113_HFC134a_fraction(:,end-15:end),16)),mean(prctile(CFC_113_HFC134a_fraction(:,end-15:end),50)),mean(prctile(CFC_113_HFC134a_fraction(:,end-15:end),84))

CFC_113_HFC125_fraction = mean(e_HFC125_CFC113(ResampleIndex_CFC113,1:end-3)./mean(emissions_CFC113(ResampleIndex_CFC113,1:end-2)));
% 2015-2019
mean(prctile(CFC_113_HFC125_fraction(:,end-4:end),16)),mean(prctile(CFC_113_HFC125_fraction(:,end-4:end),50)),mean(prctile(CFC_113_HFC125_fraction(:,end-4:end),84))
% 2004-2019
mean(prctile(CFC_113_HFC125_fraction(:,end-15:end),16)),mean(prctile(CFC_113_HFC125_fraction(:,end-15:end),50)),mean(prctile(CFC_113_HFC125_fraction(:,end-15:end),84.2))

CFC_114_bank_fraction = e_banks_CFC114(ResampleIndex_all_joint,1:end-2)./(emissions_CFC114_total(ResampleIndex_all_joint,1:end-2));
% 2015-2019
mean(prctile(CFC_114_bank_fraction(:,end-4:end),16)),mean(prctile(CFC_114_bank_fraction(:,end-4:end),50)),mean(prctile(CFC_114_bank_fraction(:,end-4:end),84))
% 2004-2019
mean(prctile(CFC_114_bank_fraction(:,end-15:end),16)),mean(prctile(CFC_114_bank_fraction(:,end-15:end),50)),mean(prctile(CFC_114_bank_fraction(:,end-15:end),84))

CFC_114_HFC134a_fraction = (e_FS_CFC114(ResampleIndex_all_joint,1:end-3)./(emissions_CFC114_total(ResampleIndex_all_joint,1:end)));
% 2015-2019
mean(prctile(CFC_114_HFC134a_fraction(:,end-4:end),16)),mean(prctile(CFC_114_HFC134a_fraction(:,end-4:end),50)),mean(prctile(CFC_114_HFC134a_fraction(:,end-4:end),84))
% 2004-2019
mean(prctile(CFC_114_HFC134a_fraction(:,end-15:end),16)),mean(prctile(CFC_114_HFC134a_fraction(:,end-15:end),50)),mean(prctile(CFC_114_HFC134a_fraction(:,end-15:end),84))

CFC_114_HFC125_fraction = mean(e_HFC125_CFC114(ResampleIndex_all_joint,1:end-3)./mean(emissions_CFC114_total(ResampleIndex_all_joint,1:end-2)));
% 2004-2019
mean(prctile(CFC_114_HFC125_fraction(:,end-15:end),16)),mean(prctile(CFC_114_HFC125_fraction(:,end-15:end),50)),mean(prctile(CFC_114_HFC125_fraction(:,end-15:end),84))

CFC_115_HFC125_fraction = (e_HFC125_CFC115(ResampleIndex_CFC115,1:end-3)./(emissions_CFC115(ResampleIndex_CFC115,1:end-2)));
% 2015-2019
mean(prctile(CFC_115_HFC125_fraction(:,end-4:end),16)),mean(prctile(CFC_115_HFC125_fraction(:,end-4:end),50)),mean(prctile(CFC_115_HFC125_fraction(:,end-4:end),84))
% 2004-2019
mean(prctile(CFC_115_HFC125_fraction(:,end-15:end),16)),mean(prctile(CFC_115_HFC125_fraction(:,end-15:end),50)),mean(prctile(CFC_115_HFC125_fraction(:,end-15:end),84))

HCFC_133a_HFC134a_fraction = mean(e_HFC134a_HCFC133a(ResampleIndex_all_joint,1:end-3)./mean(emissions_HCFC133a_total(ResampleIndex_all_joint,1:end-3)));
% 2015-2019
mean(prctile(HCFC_133a_HFC134a_fraction(:,end-15:end),16)),mean(prctile(HCFC_133a_HFC134a_fraction(:,end-15:end),50)),mean(prctile(HCFC_133a_HFC134a_fraction(:,end-15:end),84))
%% pathway fraction
% already resampled
pathway_1_fraction_global = pathway_1_HFC134a./(pathway_1_HFC134a+pathway_2_HFC134a);
1-mean(prctile(pathway_1_fraction_global(:,15:30),16)),1-mean(prctile(pathway_1_fraction_global(:,15:30),50)),1-mean(prctile(pathway_1_fraction_global(:,15:30),84))
1-mean(prctile(pathway_fraction_A5(ResampleIndex_all_joint,15:30),16)),1-mean(prctile(pathway_fraction_A5(ResampleIndex_all_joint,15:30),50)),1-mean(prctile(pathway_fraction_A5(ResampleIndex_all_joint,15:30),84))
1-mean(prctile(pathway_fraction_nonA5(ResampleIndex_all_joint,15:30),16)),1-mean(prctile(pathway_fraction_nonA5(ResampleIndex_all_joint,15:30),50)),1-mean(prctile(pathway_fraction_nonA5(ResampleIndex_all_joint,15:30),84))

%% fraction of CFCs produced in A 5
mean(prctile(A5_113_FS_prod(ResampleIndex_CFC113,26:30)./total_113_FS_prod(ResampleIndex_CFC113,26:30),16))
mean(prctile(A5_114_FS_prod(ResampleIndex_all_joint,26:30)./total_114_FS_prod(ResampleIndex_CFC113,26:30),16))



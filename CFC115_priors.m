function [longProdsamps,shortProdsamps,RFlongBank,RFshortBank,Years, short_fraction] = CFC115_priors(N)

%% CFC-115

load('CFC115Production.mat'); 

prod_s = CFC115Production.short; 
prod_l = CFC115Production.long; 

total1 = prod_s + prod_l; 
total2 = 0.001*CFC115Production.nonA5+0.001*CFC115Production.A5;

years1 = CFC115Production.AFEASyears; 
years2 = CFC115Production.UNEPyears;

yr1 = min(years2(2), years1(1)); 
yr2 = max(years2(end), years1(end)); 

Years = yr1:yr2;
TotalProd = zeros(length(Years),3); 

yr1 = find(Years == years1(1)); 
yr2 = find(Years == years1(end)); 
TotalProd(yr1:yr2,1) = total1; 

yr1 = find(Years == years2(2)); 
yr2 = find(Years == years2(end)); 
TotalProd(yr1:yr2,2) = total2(2:end); 

TotalProd(:,3) = max(TotalProd(:,1), TotalProd(:,2)); 

ind = Years < years1(1); 
prod_s = [zeros(sum(ind),1); prod_s];
prod_l = [zeros(sum(ind),1); prod_l];

prod_s(ind) = prod_s(sum(ind)+1); 
prod_l(ind) = prod_l(sum(ind)+1); 

ind = Years > years1(end);
prod_s = [prod_s; zeros(sum(ind),1)];
prod_l = [prod_l; zeros(sum(ind),1)];

ind2 = find(Years == years1(end));
prod_s(ind) = prod_s(ind2); 
prod_l(ind) = prod_l(ind2); 

tmp_Total = prod_s + prod_l; 
Scaling = TotalProd(:,3)./tmp_Total;

prod_s = prod_s.*Scaling; 
prod_l = prod_l.*Scaling; 

CFC115ReportedProd.BankedProd =  prod_s + prod_l;
CFC115ReportedProd.yrs = Years;
%save('Prior_Production.mat','CFC115ReportedProd','-append')

% Short Production Prior Samps
short_fraction = 0.50+0.4*rand(N,1);
long_fraction = 1 - short_fraction;
prod_short = short_fraction.*prod_l'; 
prod_long = long_fraction.*prod_l';

Nyrs = length(Years); 
rho(1,1,:) = 0.5+0.5*betarnd(2,2,N,1);
N_s = 100000; % making sampling tractable
exp_val = repmat(abs(repmat([1:Nyrs],Nyrs,1)-repmat([1:Nyrs]',1,Nyrs)),1,1,N_s);
%size(exp_val)

if N>N_s
    
    for ii = N_s:N_s:N
        rho_tmp = repmat(rho(1,1,1+(ii - N_s):ii),Nyrs,Nyrs,1);
        %ii
        %N_s
        %size(rho_tmp)
        Rhom    = rho_tmp.^exp_val;
        clear rho_tmp 
        Zm = mvnrnd(zeros(N_s,Nyrs), Rhom, N_s);
        clear Rhom
        Um = normcdf(Zm,0,1);
        clear Zm
        shortProdsamps(1+(ii - N_s):ii,1:Nyrs,:) = 0.2*prod_short(1+(ii - N_s):ii, :).*logninv(Um,zeros(N_s,Nyrs),0.5*ones(N_s,Nyrs))+0.80*prod_short(1+(ii - N_s):ii,:);
    end
end

Nyrs = length(Years); 
rho(1,1,:) = 0.5+0.5*betarnd(2,2,N,1);
exp_val = repmat(abs(repmat([1:Nyrs],Nyrs,1)-repmat([1:Nyrs]',1,Nyrs)),1,1,N_s);
if N>N_s
    
    for ii = N_s:N_s:N
        rho_tmp = repmat(rho(1, 1, 1+(ii - N_s):ii),Nyrs,Nyrs,1);
        Rhom    = rho_tmp.^exp_val;
        clear rho_tmp 
        Zm = mvnrnd(zeros(N_s,Nyrs), Rhom, N_s);
        clear Rhom
        Um = normcdf(Zm,0,1);
        clear Zm
        longProdsamps(1+(ii - N_s):ii,1:Nyrs) = (0.2*prod_long(1+(ii - N_s):ii, :)).*logninv(Um,zeros(N_s,Nyrs),0.5*ones(N_s,Nyrs))+(0.8*prod_long(1+(ii - N_s):ii, :));
    end
end


RFshortBank.y1 = betarnd(8,8,1,N); 
RFshortBank.y2 = 1 - RFshortBank.y1; 


RFlongBank.y1 = 2*0.07*betarnd(4,4,1,N); 

% update from Megan
% m = 1 - exp(-1/10); v = (m)^2;
% mu = log(m^2/sqrt(v+m^2));
% sigma = sqrt(log(v/m^2+1));
% RFlongBank.LT = lognrnd(mu,sigma,1,N);

% original used in her 2022 paper
m = 1 - exp(-1/10); v = (0.50*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
RFlongBank.LT = lognrnd(mu,sigma,1,N);


%save('CFC115Priors.mat','longProdsamps','shortProdsamps','RFlongBank','RFshortBank','Years');



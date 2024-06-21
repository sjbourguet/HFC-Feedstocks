function [longProdsamps, shortProdsamps, RFlongBank, RFshortBank, Years] = CFC114_priors(N)
%% CFC-114
load('CFC114Production.mat'); 

prod_s = CFC114Production.short; 
prod_l = CFC114Production.long; 

total1 = prod_s + prod_l; 
total2 = 0.001*CFC114Production.nonA5+0.001*CFC114Production.A5;

years1 = CFC114Production.AFEASyears; 
years2 = CFC114Production.UNEPyears;

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
%size(prod_s)
CFC114ReportedProd.BankedProd =  prod_s + prod_l;
CFC114ReportedProd.yrs = Years;
%save('Prior_Production.mat','CFC114ReportedProd','-append')

% Short Production Prior Samps
Nyrs = length(Years); 
rho(1,1,:) = 0.5+0.5*betarnd(2,2,N,1);
N_s = 100000; % making sampling tractable
exp_val = repmat(abs(repmat([1:Nyrs],Nyrs,1)-repmat([1:Nyrs]',1,Nyrs)),1,1,N_s);

if N>N_s
    
    for ii = N_s:N_s:N
        rho_tmp = repmat(rho(1,1,1+(ii - N_s):ii),Nyrs,Nyrs,1);
        Rhom    = rho_tmp.^exp_val;
        clear rho_tmp 
        Zm = mvnrnd(zeros(N_s,Nyrs), Rhom, N_s);
        clear Rhom
        Um = normcdf(Zm,0,1);
        clear Zm
        shortProdsamps(1+(ii - N_s):ii,1:Nyrs) = 0.2*repmat(prod_s',N_s,1).*logninv(Um,zeros(N_s,Nyrs),0.5*ones(N_s,Nyrs))+0.95*repmat(prod_s',N_s,1);
    end
end


rho(1,1,:) = 0.5+0.5*betarnd(2,2,N,1);
N_s = 100000; % making sampling tractable
exp_val = repmat(abs(repmat([1:Nyrs],Nyrs,1)-repmat([1:Nyrs]',1,Nyrs)),1,1,N_s);

if N>N_s
    
    for ii = N_s:N_s:N
        rho_tmp = repmat(rho(1,1,1+(ii - N_s):ii),Nyrs,Nyrs,1);
        Rhom    = rho_tmp.^exp_val;
        clear rho_tmp 
        Zm = mvnrnd(zeros(N_s,Nyrs), Rhom, N_s);
        clear Rhom
        Um = normcdf(Zm,0,1);
        clear Zm
        longProdsamps(1+(ii - N_s):ii,1:Nyrs) = 0.2*repmat(prod_l',N_s,1).*logninv(Um,zeros(N_s,Nyrs),0.5*ones(N_s,Nyrs))+0.95*repmat(prod_l',N_s,1);
    end
end

%longProdsamps(:,1:Nyrs) = (repmat(0.2*prod_l',N,1)).*logninv(Um,zeros(N,Nyrs),0.5*ones(N,Nyrs))+0.90*repmat(prod_l',N,1);

RFshortBank.y1 = 0.4 + 0.2*betarnd(8,8,1,N); 
RFshortBank.y2 = 1 - RFshortBank.y1; 

m = 0.02; v = (m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
RFlongBank.y1 = lognrnd(mu,sigma,1,N);

m = 1 - exp(-1/20); v = (0.50*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
RFlongBank.LT = lognrnd(mu,sigma,1,N);

%save('CFC114Priors.mat','longProdsamps','shortProdsamps','RFlongBank','RFshortBank','Years');


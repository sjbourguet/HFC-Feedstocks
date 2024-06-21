function [b,e] = Bank_Emiss_Conc(ProdIn,RFIn,DEIn,BanksIn,LT)
if size(DEIn,1) == 1
    DEIn = DEIn';
end
if size(RFIn,1) == 1
    RFIn = RFIn';
end

    b = ProdIn-DEIn.*ProdIn+(ones(size(RFIn))-RFIn).*BanksIn;
    e = RFIn.*BanksIn+DEIn.*ProdIn;

end
function [pX] = BPE_uncertainty_plots(Var, Year, color_var)

LB = prctile(Var,18);
MED = prctile(Var,50);
UB  = prctile(Var,84);
pX = boundedline(Year, MED',[MED'-LB',UB'-MED'],'alpha','cmap',color_var);
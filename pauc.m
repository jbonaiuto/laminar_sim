function [pvalue Wold Wnew]= pauc(predOld,predNew,outcome)
% Compares two AUCs derived from same cases
% Instead author recommands the use of:
%       - NetReclassificationImprovement.m
%       - IntegratedDiscriminationImprovement.m
%
% Usage:
%    [pvalue Wold Wnew] = pauc(predOld,predNew,outcome)
%
%
% (c) Louis Mayaud, 2011 (louis.mayaud@gmail.com) 
% Please reference :
% Mayaud, Louis, et al. "Dynamic Data During Hypotensive Episode Improves
%        Mortality Predictions Among Patients With Sepsis and Hypotension*."
%        Critical care medicine 41.4 (2013): 954-962.

[Wold SEold] = wilcoxon_satistic(predOld,outcome);
[Wnew SEnew] = wilcoxon_satistic(predNew,outcome);

Revent = abs(corr(predOld(outcome==1),predNew(outcome==1),'type','Kendall'));
RNevent = abs(corr(predOld(outcome==0),predNew(outcome==0),'type','Kendall'));

Ravg = (Revent + RNevent)/2 ;
Wavg = (Wold + Wnew)/2;

load Rcoeff-AUC-sign.mat
AUCs = .7:.025:.975;
Rs = .02:.02:.9;

[m RIdx] = min(sqrt((Rs-Ravg).^2));
[m WIdx] = min(sqrt((AUCs-Wavg).^2));
R = correlationAUC(RIdx,WIdx); % Get the value from the table 
    % See Table I in Hanley et Mc Neil, A method of Comparing the AUROC derived from same cases
    % Radiology vol. 148, N. 3, pp839-843, Sept. 1983
    
    if SEold^2 + SEnew^2 - 2*R*SEold*SEnew <0
        pvalue = NaN;
        return
    end
    
z = abs(Wnew-Wold)/sqrt( SEold^2 + SEnew^2 + 2*R*SEold*SEnew ) ;
pvalue = Ztest(z);

function [W SE] = wilcoxon_satistic(pred,outcome)
% from Hanley et McNeil Apr. 1982
% Radiology 143, The meaning and use of the area under the ROC curve
% Returns Wilcoxon Stat (equivelent to AUC) and its standard error

X = pred(outcome==1);
Y = pred(outcome==0);

x=length(X);
y = length(Y);


Xa = pred(outcome==1);
Xn = pred(outcome==0);
w=0;
q1=0; q2 = 0;
for i=1:length(Xa)
    for j=1:length(Xn)
        if Xa(i) > Xn(j)
            w = w+1;
        elseif Xa(i)==Xn(j)
            w = w + 1/2;
        end
        
    end
end
W = w/(length(Xa)*length(Xn)) ; 
q1=W/(2-W);
q2=2*W^2/(1+W);
SE = sqrt( (W*(1-W)+(x-1)*(q1-W^2) + (y-1)*(q2-W^2)) / (x*y) ) ;   

function p=Ztest(Z)
% From statbag pnorm.m (Steve Roberts)
% Return probability of z > Z if z is distributed as N(0,1)
root2=sqrt(2);
x=-1*Z/root2;
p=1-0.5*erfc(x);
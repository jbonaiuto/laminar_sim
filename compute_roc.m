function [fpr,tpr,thresholds,auc]=compute_roc(scores, labels, Nsims)
    scores(end+1)=0.0;
    [sScores,sorted] = sort(scores,'descend');
    thresholds=sScores;    
    tp=zeros(1,length(thresholds));
    fp=zeros(1,length(thresholds));
    
    for i=1:length(thresholds)
        threshold=thresholds(i);
        for idx=1:Nsims*2
            if scores(idx)>=threshold
                if labels(idx)>0
                    tp(i)=tp(i)+1;
                else
                    fp(i)=fp(i)+1;
                end
            end
        end
    end
    tpr=tp./(Nsims);
    fpr=fp./(Nsims);
    auc=trapz(fpr,tpr);
end
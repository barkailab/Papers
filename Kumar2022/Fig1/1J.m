%% prom Idx without subtel genes

% needs whichProm, promVecIdx, GP, CISBP, medianMotif, scgenome, nmerRed
subtelomereGenes = GP.groups{1, 23}{1, 2}{1, 61};
promoterIDXvec(ismember(whichProm,subtelomereGenes))=0;
descTF = readtable('/home/labs/barkailab/divyakr/matlab/FinalForThesis/descTF.xlsx');
myCurrTable = descTF;
myCurrTable(cellfun(@isempty,myCurrTable.DBDname),:)=[];
myCurrTable(ismember(myCurrTable.FAMILY, 'Zinc cluster'),:)= [];
myCurrTable(ismember(myCurrTable.TF, 'HOT1'),:)= [];
myCurrTable(ismember(myCurrTable.TF, 'MSN1'),:)= [];
 myCurrTable(ismember(myCurrTable.TF, {'PIP2';'ARO80';'CAT8';'CHA4';'PHD1';'MSN1';'GCR1'}),:)= [];
myCurrTable.DBDname(ismember(myCurrTable.TF, 'YAP1'))= {'YAP1_d2_500'};
for TFNo=1:length(myCurrTable.TF)
    TFReq = myCurrTable.TF{TFNo}; DBDName = myCurrTable.DBDname{TFNo};
    
    [ ~, topIVMIdxall] = maxk(CISBP7merDist.(TFReq), 10);
    [ ~, topIVMIdxWT] = max(medianMotifNew.(TFReq)(topIVMIdxall));  topIVMIdxWT = topIVMIdxall( topIVMIdxWT );
    [ ~, topIVMIdxDBD] = max(medianMotifNew.(DBDName)(topIVMIdxall));  topIVMIdxDBD = topIVMIdxall( topIVMIdxDBD ); clearvars topIVMIdxall;
    
    [a,b] =regexp(SC_genome, nmerRed7.seq(topIVMIdxWT), 'start', 'end');
    [c,d] =regexp(SC_genome, nmerRed7.rcSeq(topIVMIdxWT), 'start', 'end');
    motifPosWT = intersect(find(promoterIDXvec), round(mean([[cat(2 , a{:}, c{:})] ; [cat(2 , b{:}, d{:})]])));
    clearvars a b c d topIVMIdxWT
    
    [a,b] =regexp(SC_genome, nmerRed7.seq(topIVMIdxDBD), 'start', 'end');
    [c,d] =regexp(SC_genome, nmerRed7.rcSeq(topIVMIdxDBD), 'start', 'end');
    motifPosDBD = intersect(find(promoterIDXvec), round(mean([[cat(2 , a{:}, c{:})] ; [cat(2 , b{:}, d{:})]])));
    clearvars a b c d topIVMIdxDBD
    motifPosWT = [motifPosWT motifPosDBD]; motifPosDBD =motifPosWT;
    
    WT = load(['/home/labs/barkailab/divyakr/matlab/FinalForThesis/newFiles/structs/medianNorm/' TFReq '.mat']);
    DBD = load(['/home/labs/barkailab/divyakr/matlab/FinalForThesis/newFiles/structs/medianNorm/' DBDName '.mat']);
    MNase = load('/home/labs/barkailab/divyakr/matlab/FinalForThesis/newFiles/structs/medianNorm/TDH3P_MNase.mat');
    rescaledSmWT = movmean((WT.medianNorm./sum(WT.medianNorm(promoterIDXvec==1)))*sum(promoterIDXvec==1), 51);
    rescaledSmDBD = movmean((DBD.medianNorm./sum(DBD.medianNorm(promoterIDXvec==1)))*sum(promoterIDXvec==1), 51);
    rescaledSmMNase = movmean((MNase.medianNorm./sum(MNase.medianNorm(promoterIDXvec==1)))*sum(promoterIDXvec==1), 51);
    enr.name(TFNo) = {TFReq};
    enr.WT(TFNo) = mean(rescaledSmWT(motifPosWT))/mean(rescaledSmWT(promoterIDXvec==1));
    enr.DBD(TFNo) = mean(rescaledSmDBD(motifPosWT))/mean(rescaledSmDBD(promoterIDXvec==1));
    enr.MNase(TFNo) = mean(rescaledSmMNase(motifPosWT))/mean(rescaledSmMNase(promoterIDXvec==1));
    clearvars WT DBD rescaledSmDBD rescaledSmWT motifPosDBD motifPosWT TFReq DBDName MNase rescaledSmMNase
end

scatter( log2(enr.WT), log2(enr.DBD),200,'filled'); xlim(ylim);%xlim([0 70]) , 'MarkerFaceColor', [0.0069    0.3040    0.5961]
 hold on; plot([log2(2) log2(2)], ylim, '--k');
 hold on; plot(xlim, [log2(2) log2(2)], '--k')
set(gcf, 'Renderer', 'painters', 'Color', [1 1 1]);
xlabel('log2(motif  Enr WT)'); ylabel('log2(motifEnr DBD)')

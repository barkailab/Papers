descTF = readtable('/home/labs/barkailab/divyakr/matlab/FinalForThesis/descTF.xlsx');
myCurrTable = descTF;
myCurrTable(cellfun(@isempty,myCurrTable.DBDname),:)=[];
myCurrTable(ismember(myCurrTable.FAMILY, 'Zinc cluster'),:)= [];
myCurrTable(ismember(myCurrTable.TF, 'HOT1'),:)= [];
myCurrTable(ismember(myCurrTable.TF, 'MSN1'),:)= [];
 myCurrTable(ismember(myCurrTable.TF, {'PIP2';'ARO80';'CAT8';'CHA4';'PHD1';'MSN1';'GCR1'}),:)= [];
 myCurrTable.DBDname(ismember(myCurrTable.TF, 'YAP1'))= {'YAP1_d2_500'};

promoterIDXvec(ismember(whichProm,subtelomereGenes))=0;
for TFNo = 1:length(myCurrTable.TF)
    
    TFReq = myCurrTable.TF{TFNo}; DBDName = myCurrTable.DBDname{TFNo};
    
      [ ~, topIVMIdxall] = maxk(CISBP7merDist.(TFReq), 10);
    [~, topIVMIdxWT] = max(medianMotifNew.(TFReq)(topIVMIdxall)); topIVMIdxWT = topIVMIdxall( topIVMIdxWT );
    [~, topIVMIdxDBD] = max(medianMotifNew.(DBDName)(topIVMIdxall)); topIVMIdxDBD = topIVMIdxall( topIVMIdxDBD ); clearvars topIVMIdxall
   
    [a,b] =regexp(SC_genome, nmerRed7.seq(topIVMIdxWT), 'start', 'end');
    [c,d] =regexp(SC_genome, nmerRed7.rcSeq(topIVMIdxWT), 'start', 'end');
    motifPosWT = intersect(find(promoterIDXvec), round(mean([[cat(2 , a{:}, c{:})] ; [cat(2 , b{:}, d{:})]]))); 
    currWTlength = length(motifPosWT);
     clearvars a b c d topIVMIdxWT
    
    [a,b] =regexp(SC_genome, nmerRed7.seq(topIVMIdxDBD), 'start', 'end');
    [c,d] =regexp(SC_genome, nmerRed7.rcSeq(topIVMIdxDBD), 'start', 'end');
    motifPosDBD = intersect(find(promoterIDXvec), round(mean([[cat(2 , a{:}, c{:})] ; [cat(2 , b{:}, d{:})]]))); 
    currDBDlength = length(motifPosDBD);
    clearvars a b c d topIVMIdxDBD
    
    motifPosWT = [motifPosWT motifPosDBD]; motifPosDBD =motifPosWT;
     posMatWT = motifPosWT'+  (-25:25); posMatDBD =posMatWT;
     
    %uniformity score 1
    WT = load(['/home/labs/barkailab/divyakr/matlab/FinalForThesis/newFiles/structs/medianNorm/' TFReq '.mat']);
    DBD = load(['/home/labs/barkailab/divyakr/matlab/FinalForThesis/newFiles/structs/medianNorm/' DBDName '.mat']);
    WT.loci = sum(WT.medianNorm(posMatWT),2);  DBD.loci = sum(DBD.medianNorm(posMatDBD),2);
    UScore.name(TFNo,1) = {TFReq};
    UScore.maxLocusWT(TFNo,1) =  max( WT.loci);
     UScore.maxLocusDBD(TFNo,1) =  max( DBD.loci);
    UScore.WT(TFNo,1) =  (sum(WT.loci >(max( WT.loci)*0.03)) /length(WT.loci))*100;
    UScore.DBD(TFNo,1) =  (sum(DBD.loci >(max( DBD.loci)*0.03)) /length(DBD.loci))*100;
    UScore.noOfLociWT(TFNo,1) = length(motifPosWT); 
    UScore.noOfLociDBD(TFNo,1) = length(motifPosDBD); clearvars WT DBD;
    
    UScore.nDBD(TFNo,1) =currDBDlength;
    UScore.nWT(TFNo,1) =currWTlength;
    clearvars  topIVMIdx TFReq DBDName posMatWT posMatDBD motifPosWT motifPosDBD currWTlength currDBDlength
end
% UScore = struct2table(UScore);
%%

highlightName ={'SKO1';'GLN3'};
highlightIdx = find(contains(UScore.name, highlightName));

figure( 'Color', [1 1 1], 'Renderer', 'painters' )%'Position', [932 537 449 347],
% axes('Position', [0.15 0.15 0.75 0.75])
scatter(UScore.WT, UScore.DBD, 200, log2(UScore.noOfLociDBD),  'filled')
% addMrkr(UScore.WT', UScore.DBD', 0.015 , (UScore.noOfLociDBD),  (UScore.noOfLociWT))
text(UScore.WT(highlightIdx)+1.3, UScore.DBD(highlightIdx)+1.3, UScore.name(highlightIdx),  'FontWeight' , 'bold');
colormap(flipud(brewermap(100, 'PuBu')));cb= colorbar; ylabel(cb, 'log2(no of Loci)'); %caxis([1 max([log10(UScore.noOfLociWT);log10(UScore.noOfLociDBD)] )])
xlabel('WT'); ylabel('DBD'); %title('fraction loci> 3%max')
ylim([min([ xlim ylim]) max([ xlim ylim])]); xlim(ylim);
hold on; plot(xlim, ylim, '--k')
% caxis([0 1000]); ylabel(cb, 'number of Loci')
set(gcf, 'Renderer', 'painters', 'Color', [1 1 1], 'Position', [680 558 560 420])%, 'Position'  , [2185 292 828 672]);
ylim([-2 70]); %ylim([-0.02 0.7])
xlim(ylim); caxis([7 13])
%%
clearvars cb UScore
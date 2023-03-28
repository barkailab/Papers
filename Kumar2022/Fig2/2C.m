%%
promoterIDXvec(ismember(whichProm,subtelomereGenes))=0;
TFReq = 'SKO1'; 
DBDName = char(myCurrTable.DBDname(ismember(myCurrTable.TF, TFReq)));

      [ ~, topIVMIdxall] = maxk(CISBP7merDist.(TFReq), 10);
    [~, topIVMIdxWT] = max(medianMotifNew.(TFReq)(topIVMIdxall)); topIVMIdxWT = topIVMIdxall( topIVMIdxWT );
    [~, topIVMIdxDBD] = max(medianMotifNew.(DBDName)(topIVMIdxall)); topIVMIdxDBD = topIVMIdxall( topIVMIdxDBD ); clearvars topIVMIdxall
   
    [a,b] =regexp(SC_genome, nmerRed7.seq(topIVMIdxWT), 'start', 'end');
    [c,d] =regexp(SC_genome, nmerRed7.rcSeq(topIVMIdxWT), 'start', 'end');
    motifPosWT = intersect(find(promoterIDXvec), round(mean([[cat(2 , a{:}, c{:})] ; [cat(2 , b{:}, d{:})]]))); 
    posMatWT = motifPosWT'+  (-25:25); clearvars a b c d topIVMIdxWT
    
    [a,b] =regexp(SC_genome, nmerRed7.seq(topIVMIdxDBD), 'start', 'end');
    [c,d] =regexp(SC_genome, nmerRed7.rcSeq(topIVMIdxDBD), 'start', 'end');
    motifPosDBD = intersect(find(promoterIDXvec), round(mean([[cat(2 , a{:}, c{:})] ; [cat(2 , b{:}, d{:})]]))); 
    posMatDBD = motifPosDBD'+  (-25:25); clearvars a b c d topIVMIdxDBD
    
    %uniformity score 1
    WT = load(['/home/labs/barkailab/divyakr/matlab/FinalForThesis/newFiles/structs/medianNorm/' TFReq '.mat']);
    DBD = load(['/home/labs/barkailab/divyakr/matlab/FinalForThesis/newFiles/structs/medianNorm/' DBDName '.mat']);
    WTloci = sum(WT.medianNorm(posMatWT),2);  DBDloci = sum(DBD.medianNorm(posMatDBD),2);
    clearvars  WT DBD
     figure; 
    axes('position', [0.1300 0.575 0.7750 0.35])
    histogram(log2(WTloci), 'FaceColor', [0    0.3941    0.2100], 'BinWidth', 0.25 ); %hold on; plot([log2(0.03*max(WTloci)) log2(0.03*max(WTloci))], ylim, ':k', 'LineWidth', 5);
    hold on; fill([log2(0.03*max(WTloci)) max(xlim) max(xlim) log2(0.03*max(WTloci))], [0 0 max(ylim) max(ylim)], 'r', 'FaceAlpha', 0.5);
ylabel('no of loci'); title([TFReq]); xticks([]) %    xlabel('log2(locus signal)');
    axes('position', [0.1300 0.2 0.7750 0.35])
    histogram(log2(DBDloci), 'FaceColor', [0.7400    0.8938    0.5854], 'BinWidth', 0.25); %hold on; plot([log2(0.03*max(DBDloci)) log2(0.03*max(DBDloci))], ylim, ':k', 'LineWidth', 5);
 hold on; fill([log2(0.03*max(DBDloci)) max(xlim) max(xlim) log2(0.03*max(DBDloci))], [0 0 max(ylim) max(ylim)], 'r', 'FaceAlpha', 0.5);
    ylabel('no of loci');  xlabel('locus binding'); %title([TFReq ' DBD']);   
    set(gcf, 'Renderer', 'painters', 'Color', [1 1 1], 'Position', [785 650 405 390]); linkaxes
    clearvars DBDloci DBDName motifPosDBD motifPosWT posMatDBD posMatWT WTloci
descTF = readtable('/home/labs/barkailab/divyakr/matlab/FinalForThesis/descTF.xlsx');
myCurrTable = descTF;
myCurrTable(cellfun(@isempty,myCurrTable.DBDname),:)=[];

%required: nucleosome  median, CISBP, nmerRed7, GP, Scgenome

TFReq = 'ACE2'; geneTarget = 'YLR049C';
 TFNo = find(ismember(myCurrTable.TF, TFReq)); currStrains = {TFReq; myCurrTable.DBDname{TFNo};  myCurrTable.nonDBDName{TFNo}};
 h3smooth = NucleosomeData.median';
 
[ ~, topIVMIdxall] = maxk(CISBP7merDist.(TFReq), 10);
% [~, x] = max(medianMotifNew.(TFReq)(topIVMIdxall)); topIVMIdx(1) = topIVMIdxall( x ); clearvars x %new way
[~, x] = max(medianMotifNew.(currStrains{2})(topIVMIdxall)); topIVMIdx(1) = topIVMIdxall( x );
% topIVMIdx = topIVMIdxall;
useMotif = [nmerRed7.seq{(topIVMIdx)} '|' nmerRed7.rcSeq{(topIVMIdx)}]; clearvars x %new way


load('/home/labs/barkailab/divyakr/matlab/FinalForThesis/medianTemps/promoterLengthsORF.mat');
load('/home/labs/barkailab/divyakr/matlab/FinalForThesis/medianFinals/promoterIDXvec.mat')

g=ismember(GP.gene_infoR64.name, geneTarget);
posG= GP.gene_infoR64.position(g,:);
promoterPos = GP.chrIdx(posG(1))+ posG(2)+ ((-promoterLengthsORF(g):100).*GP.gene_infoR64.dir(g)) ;

load(['/home/labs/barkailab/divyakr/matlab/FinalForThesis/newFiles/structs/medianNorm/'  currStrains{1} '.mat'])
imageMat(1,:) = medianNorm(promoterPos,:)';
maxVal(1) = max(medianNorm(logical(promoterIDXvec)));
clearvars medianNorm
load(['/home/labs/barkailab/divyakr/matlab/FinalForThesis/newFiles/structs/medianNorm/'  currStrains{2} '.mat'])
imageMat(2,:) = medianNorm(promoterPos,:)';
maxVal(2) = max(medianNorm(logical(promoterIDXvec)));
clearvars medianNorm
load(['/home/labs/barkailab/divyakr/matlab/FinalForThesis/newFiles/structs/medianNorm/'  currStrains{3} '.mat'])
imageMat(3,:) = medianNorm(promoterPos,:)';
maxVal(3) = max(medianNorm(logical(promoterIDXvec)));
clearvars medianNorm


%adjusting promoter seq accoring to direction
if GP.gene_infoR64.dir(g)>0
    promoterSeq = SC_genome(promoterPos);
else
    promoterSeq = seqrcomplement(SC_genome(sort(promoterPos)));
end
promoterSeq = upper(promoterSeq);

% %find motif positions
motifStr = regexp(char(useMotif), '[|]', 'split');
[motStart1, motEnd1] = regexp(promoterSeq, char(motifStr{1}), 'start', 'end');
[motStart2, motEnd2] = regexp(promoterSeq, char(motifStr{2}), 'start', 'end');
motMid = mean([[motStart1 motStart2]; [motEnd1 motEnd2]])-promoterLengthsORF(g);
%find motif positions
% motifStr =useMotif;
% for i=1:7
%    if ~isempty( regexp(promoterSeq, char(useMotif{i})))
% [motStart1(i), motEnd1(i)] = regexp(promoterSeq, char(useMotif{i}), 'start', 'end');
%    end
% end
% motMid = mean([[motStart1]; [motEnd1]])-promoterLengthsORF(g);


figure('Position', [978 239 492 603], 'Color', [1 1 1], 'Renderer', 'painters')
for i = 1:3
    axes('Position', [0.22 0.4500-(0.11*(i-1)) 0.7 0.11])
    area((-promoterLengthsORF(g):100), rescale(h3smooth(promoterPos,:), 0, 1, 'InputMin', 0, 'InputMax', 300), 'FaceColor', [0.93 0.93 0.88], 'EdgeColor', 'none');
%     title(strrep(currStrains(i), '_', ' '), 'Position', [-promoterLengthsORF(g) 0.35 0], 'HorizontalAlignment', 'right' )
    hold on
    plot((-promoterLengthsORF(g):100), rescale(imageMat(i,:), 0, 1, 'InputMin', 0, 'InputMax', 0.6*maxVal(i)) , 'Color', [0 0 0]);
    hold on; sc = scatter(motMid, zeros(size(motMid))+0.03, 50,  'filled', 'o', 'DisplayName', motifStr{1}, 'MarkerFaceColor',  [0.7096    0.0133    0.0048]);
    hold on; plot([0 0], [0 1], 'k:')
    yticks([])
    title(strrep(currStrains(i), '_', ' '), 'Position', [-701 0.35 0], 'HorizontalAlignment', 'right' )
    xlim([-700 100]); ylim([0 1]) %xlim([-promoterLengthsORF(g) 100]), ylim([0 1])
    if i==1
%         text(100-promoterLengthsORF(g)/2, 1.1, GP.gene_infoR64.name(g), 'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')
         text(100-700/2, 1.1, GP.gene_infoR64.name(g), 'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlignment', 'center')
    end
    if i==1||i==2
     xticks([])
    end
end
legend(sc, 'FontSize', 8, 'Position'  ,  [0.5760 0.5160 0.5000 0.1100], 'Box', 'off') ;
clearvars motStart1 motEnd1 motStart2 motEnd2 motMid


clearvars currStrains dbdEnd dbdStart g geneTarget i imageMat k maxVal motifStr posG promoterPos promoterSeq sc TFNo TFReq tops topsX topsY useMotif 

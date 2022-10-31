%% Figure 5—figure supplement 1
clearvars -except checWTdelLactisSwap
load('summaryTable.mat')
geneName = reshape([summaryTable.p1';summaryTable.p2'], 1,[]);
h3cc=load('./H3CC_henikoff.mat');
h3smooth=conv(mean(h3cc.p_all,2),normpdf([-50:50]',0,25),'same');
load('promoterLengthsORF.mat')
load('promoterIDXvec.mat');
load('SC_genome.mat')
GP=load('./group_imp.mat')
load('./promCorrSort.mat')
load('./paraSeqs.mat');

%% Figure 5—figure supplement 1B-E
clearvars -except GP checWTdelLactisSwap
allSamples=fieldnames(checWTdelLactisSwap.sumProm)
zscoreTH = 3.5;
quantileTH = 0.99;
Xlimit = [-4 4];
gb = true(sum(GP.chr_len(1:16)), 1);
gb(GP.chrIdx(12)+[451575-51:468958+51]) = false;
yLimitSignal = 0.5;
intBases = intBasesVec(1:6701);
intBasesMotifs = createIntBasesForMotif();
intoGene = 150;
promoterL = repmat(700,6701,1);

lineCol = lines(2);
%col = cbrewer2('Dark2',1);
col = brewermap(1,'Dark2');

cMapCorr = flipud(bone);
%cMapScatter = cbrewer2('OrRd');
cMapScatter = brewermap(256,'OrRd');

cMapScatter = cMapScatter(50:end,:);
markerCol = cMapScatter(128,:);
colSignal = repmat([0 0 0; 0.5 0.5 0.5], 2,1);
colMap = brighten(flipud(bone),0.2);

yspacer = 0.05;
xspacer = 0.05;
Wsignal = 0.25;
Wscatter = 0.04;
Wcorrmat = 0.05/2;
Hscatter =Wscatter*1.8;
Hlog2 = 0.05/2;
xPos = [0.05: 3*Wscatter+3*xspacer: 1];
xPos = repmat(xPos(1:3),5,1);
yPos = [1-(Hscatter+Hlog2+yspacer):-(Hscatter+Hlog2+2*yspacer):0];
yPos = repmat(yPos',1,3);
goodPos = createIntBasesForMotif();

figure('Units','normalized','Position',[0.5000 0.0370 0.5000 0.8907], 'color','w')
c=1;
for i = [1:5, 11:15, 21:25, 6:10, 16:20, 26:30]
    if i== 6
        save_gf(gcf,sprintf('FigS3Part%d',i),'paper','tamar21','type','svg','size',[])
        close(gcf)
        figure('Units','normalized','Position',[0.5000 0.0370 0.5000 0.8907], 'color','w')
        c=1;
    end
    selTargets = summaryTable.nHighTargetsIdx{i};
    if ~all(isnan(summaryTable.WTmutantCorr(i,:)))
        if any(isnan(summaryTable.WTmutantCorr(i,:)))
            idxChangedPar = find(~isnan(summaryTable.WTmutantCorr(i,:)));
        else
            [~, idxChangedPar] = min(summaryTable.WTmutantCorr(i,:));
        end
        F1 = summaryTable.(sprintf('p%d', idxChangedPar)){i};
        F2 = summaryTable.(sprintf('p%d', 3-idxChangedPar)){i};
        m1 = [F1, '_d', upper(F2)];
        m2 = [F2, '_d', upper(F1)];
        nRepDel = [0,0];
        currTFs = {F1,F2};
        
        clear logChange significantGenes zScoreTF
        for tf = 1:2
            currF1 = currTFs{tf};
            currF2 = currTFs{3-tf};
            currM1 = [currF1,'_d',upper(currF2)];
            zscore_WT1 = nanZscore(checWTdelLactisSwap.sumProm.(currF1));
            fitGenes =  zscore_WT1 > min(zscoreTH, quantile(zscore_WT1, quantileTH));
            if any(contains(allSamples, currM1))
                bestFit = robustfit(checWTdelLactisSwap.sumProm.(currF1)(fitGenes), checWTdelLactisSwap.sumProm.(currM1)(fitGenes), [], [],'off' );
                normFactor(tf) = log2(bestFit);
                logChange(tf,:) = log2(checWTdelLactisSwap.sumProm.(currM1)+700)- log2(checWTdelLactisSwap.sumProm.(currF1)+700)- normFactor(tf);
                significantGenes(tf,:) =  abs((logChange(tf,:)-median(logChange(tf,fitGenes)))/std(logChange(tf,fitGenes))) >= 1;
                [~, ~,~ , repIdx] = getRepeatsCorr({currM1},'dataType','sumProm');
                nRepDel(1,tf) = numel(repIdx);
            end
        end
        
        % scatter paralog that changed
        if any(contains(allSamples, m1))
            axes('Position', [xPos(c) yPos(c)+Hlog2 Wscatter Hscatter])
            maxX = max(checWTdelLactisSwap.sumProm.(F1));
            maxY = max(checWTdelLactisSwap.sumProm.(m1));
            maxCol = max(checWTdelLactisSwap.sumProm.(F2));
            scatter(checWTdelLactisSwap.sumProm.(F1)/maxX,...
                checWTdelLactisSwap.sumProm.(m1)/maxY,15,...
                checWTdelLactisSwap.sumProm.(F2)/maxCol, 'filled')
            hold on
            plot(xlim, xlim*(2^normFactor(1))*maxX/maxY, 'k:')
            xlabel([F1, '(P1)'], 'fontSize',8)
            ylabel([F1, 'Delta', upper(F2),'(m1)'], 'fontSize',8) %num2str(nRepDel(1))
            % text(max(xlim)*0.05,max(yLim)*0.965,sprintf('r = %.2f', corr(checWTdelLactisSwap.sumProm.(F1),checWTdelLactisSwap.sumProm.(m1), 'rows','pairwise')), 'fontSize',10)
            title(sprintf('r = %.2f', corr(checWTdelLactisSwap.sumProm.(F1),checWTdelLactisSwap.sumProm.(m1), 'rows','pairwise')), 'fontSize',8, 'FontWeight','normal')
            set(gca,'XTick',[0 1],'YTick',[0 1])
            cbr = colorbar()
            title(cbr, [F2, '(P2)'], 'fontSize',8)
            %cbr.Position = cbr.Position.*[1 1 0.5 1];
            caxis([0 0.8*max(caxis)])
            colormap(gca, cMapScatter)
            cbr.Position = [xPos(c)+Wscatter+0.005 yPos(c)+Hlog2 0.002 Hscatter]
            set(cbr, 'Ticks',[0,0.4,0.7]);
        end
        
        % scatter other paralog
        if any(contains(allSamples, m2))
            axes('Position', [xPos(c)+Wscatter+xspacer yPos(c)+Hlog2 Wscatter Hscatter])
                      maxX = max(checWTdelLactisSwap.sumProm.(F2));
            maxY = max(checWTdelLactisSwap.sumProm.(m2));
            maxCol = max(checWTdelLactisSwap.sumProm.(F1));
            scatter(checWTdelLactisSwap.sumProm.(F2)/maxX,...
                checWTdelLactisSwap.sumProm.(m2)/maxY,15,...
                checWTdelLactisSwap.sumProm.(F1)/maxCol, 'filled')
            hold on
            plot(xlim, xlim*(2^normFactor(2))*maxX/maxY, 'k:')
            xlabel([F2, '(P2)'], 'fontSize',8)
            ylabel([F2, 'Delta', upper(F1),'(m2)'], 'fontSize',8) %num2str(nRepDel(2))
            yLim = ylim;
            title(sprintf('r = %.2f', corr(checWTdelLactisSwap.sumProm.(F2),...
                checWTdelLactisSwap.sumProm.(m2), 'rows','pairwise')), 'fontSize',8, 'fontWeight','normal')
            ylim(yLim)
            set(gca,'XTick',[0 1],'YTick',[0 1])
            cbr = colorbar()
            title(cbr, [F1,'(P1)'], 'fontSize',8)
            caxis([0 0.8*max(caxis)])
            colormap(gca, cMapScatter)
            cbr.Position = [xPos(c)+2*Wscatter+xspacer+0.005 yPos(c)+Hlog2 0.002 Hscatter];
            set(cbr, 'Ticks',[0,0.4,0.7]);
        end
        
        % corr matrix
        axes('Position', [xPos(c)+2*(Wscatter+xspacer) yPos(c)+Hlog2 Wscatter Hscatter])
        % WITH REPEATS
        tRaw =  {F1, m1, F2, m2 };
        t = {currF1, strrep(m1, '_d', ' Delta'), currF2, strrep(m2, '_d', ' Delta')};
        [~, idxVec, sumPromRep, repIdx] = getRepeatsCorr(tRaw,'dataType','sumProm');
        [~, ~, mer7Rep, ~] = getRepeatsCorr(tRaw,'dataType','7mer');
        
        samplesWOrepeats = unique(idxVec(all(isnan(sumPromRep),1)));
        for z = samplesWOrepeats'
            if isfield(checWTdelLactisSwap.sumProm, tRaw{z})
                sumPromRep(:,idxVec==z) = repmat(checWTdelLactisSwap.sumProm.(tRaw{z}),1,sum(idxVec==z));
                currNorm = chromosomes2fullProfile(checWTdelLactisSwap, tRaw(z));
                currMer = mer_occupancy(currNorm,7,'intBases', goodPos,'method','else');
                mer7Rep(:,idxVec==z) = repmat(currMer.score,1,sum(idxVec==z));
            end
        end
        
        sumPromCorrMat = corr(sumPromRep,'rows','pairwise');
        mer7CorrMat = corr(mer7Rep,'rows','pairwise');
        combineMat = tril(sumPromCorrMat) + triu(mer7CorrMat,1);
        
        imagesc(combineMat);
        borders = [0; cumsum(accumarray(idxVec,1))]+0.5;
        tickPos = movmean(borders,2,'Endpoints','discard');
        hold on
        plot(repmat(borders',2,1), repmat(ylim',1, numel(borders)), 'k')
        plot(repmat(ylim',1, numel(borders)), repmat(borders',2,1), 'k')
        plot(xlim, ylim,'color', [0.7 0.7 0.7], 'LineWidth',1)
        
        caxis([0 1]);
        axis square
        set(gca,'XTick',[],'YTick',[]);
        colormap(gca, [0.7 0.7 0.7; colMap]);
        set(gca,'YTick', tickPos,  'YTickLabel', {'P1',sprintf('m1 (%d)',nRepDel(1)),'P2',sprintf('m2 (%d)',nRepDel(2))}, 'fontSize',8)
%         cbr = colorbar()
%         ylabel(cbr,'correlation', 'fontSize',10)
%         set(cbr, 'Ticks', [0.2:0.4:1])
        
        % log2 change
        axes('Position', [xPos(c) yPos(c)-Hlog2-0.01 3*Wscatter+2*xspacer Hlog2])
        wtVec(1,:) = checWTdelLactisSwap.sumProm.(F1);
        wtVec(2,:) = checWTdelLactisSwap.sumProm.(F2);
        
        for tf = 1:2
            plot([Xlimit], [tf tf], 'Color', [0.7 0.7 0.7], 'LineStyle', ':');
            currF1 = currTFs{tf};
            currF2 = currTFs{3-tf};
            m1 = [currF1, '_d', upper(currF2)];
            m2 = [currF2, '_d', upper(currF1)];
            hold on
            if any(contains(allSamples,m1))
                grayIdx = selTargets(~significantGenes(tf,selTargets));
                significantIdx = selTargets(significantGenes(tf,selTargets));
                sizeVec = rescale(checWTdelLactisSwap.sumProm.(currF1)(selTargets),10,90, 'InputMax', max(checWTdelLactisSwap.sumProm.(currF1)));
                % color rescaled over all genes
                fullColVec =  rescale(wtVec(3-tf,:),0,1, 'InputMax', 0.8*max(wtVec(3-tf,:)), 'InputMin',0);
                colVec = fullColVec(selTargets);
                % color rescaled only to selected targets
                %colVec = rescale(checWTdelLactisSwap.sumProm.(examplePairs{i,3-tf})(selTargets),0,1, 'InputMax', 0.8*max(wtVec(3-tf,:)), 'InputMin',0);
                scatter(logChange(tf,grayIdx), tf*ones(numel(grayIdx),1),...
                    sizeVec(ismember(selTargets,grayIdx)), [0.7 0.7 0.7],'filled',...
                    'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', [0.7 0.7 0.7]);
                scatter(logChange(tf, significantIdx), tf*ones(numel(significantIdx),1), sizeVec(ismember(selTargets,significantIdx)),...
                    colVec(ismember(selTargets,significantIdx)),...
                    'filled','MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', cMapScatter(end,:));
            end
        end
        %set(gca, 'YDir', 'reverse')
        ylim([0.5 2.5])
        plot([0 0],ylim, 'k--')
        colormap(gca, cMapScatter)
        caxis([0 1])
        xlabel('log2 fold change (DeltaParalog/wt)', 'fontSize',8, 'FontWeight','normal')
        set(gca, 'YTick', [1,2], 'YTickLabel', {F1,F2}, 'fontSize',8, 'XTick', [-4:2:4], 'TickLength',[0 0])
        
        %         % colorbar log2
        %         axes('Position', [xPos(i)+0.007 yPos-2*yspacer-Hsignal/2-0.05 0.05 0.01])
        %         colormap(gca, cMapScatter)
        %         cbr = colorbar('Location','south');
        %         cbr.AxisLocation =  'out';
        %         set(cbr, 'Ticks', [min(caxis), max(caxis)], 'TickLabels', {'min','max'}, 'TickLength', [0 0], 'FontSize',7)
        %         title(cbr, 'Paralogs signal', 'fontSize',10)
        %         axis off
        
        
        c=c+1;
    end
end
%save_gf(gcf,sprintf('FigS3Part%d',i),'paper','tamar21','type','svg','size',[])

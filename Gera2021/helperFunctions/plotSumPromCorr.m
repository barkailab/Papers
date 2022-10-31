function [SumPromCorrMat,pl] = plotSumPromCorr(StrainsList, ChecStruct, p, varargin)
% correlation of signals on promotors
%Strains = fieldnames(ChecStruct.norm);
ip = inputParser;
ip.addParameter('Type','Pearson');
ip.addParameter('axes',[]);
ip.addParameter('scale',{});

ip.parse(varargin{:});
CorrType = ip.Results.Type;
if p==1
if numel(ip.Results.axes) == 0
    figure;
else
    axes(ip.Results.axes)
end
end
if isempty(ip.Results.scale)
    myMat = nan(6701, length(StrainsList));
    for i=1:length(StrainsList)
        if isfield(ChecStruct.sumProm,StrainsList{i})
            myMat(:,i) = ChecStruct.sumProm.(StrainsList{i});
        end
    end
    SumPromCorrMat = corr(myMat,'rows','pairwise','Type',CorrType);
    
    if p==1;
        l = strrep(StrainsList,'_',' ');
        l = strrep(l,'_d',' \Delta')
        pl = imagesc(SumPromCorrMat);
        h = colorbar;
        ylabel(h,'Promoters signal corr','fontsize',15, 'FontWeight','bold');
        axis square;
        set(gca,'xtick',[1:length(StrainsList)],'XTickLabel',l,...
            'ytick',[1:length(StrainsList)],'yTickLabel',l,'fontsize',9.5);
        xtickangle(90);
        set(gcf,'color','w');
    end
else
    for i=1:length(StrainsList)
        myMat(:,i) = log2(ChecStruct.sumProm.(StrainsList{i})+700);
    end
    SumPromCorrMat = corr(myMat,'rows','pairwise','Type',CorrType);
    
    if p==1;
        l = strrep(StrainsList,'_',' ');
        l = strrep(l,'_d',' \Delta')
        pl = imagesc(SumPromCorrMat);
        h = colorbar;
        ylabel(h,'log2 Promoters corr +700','fontsize',15, 'FontWeight','bold');
        axis square;
        set(gca,'xtick',[1:length(StrainsList)],'XTickLabel',l,...
            'ytick',[1:length(StrainsList)],'yTickLabel',l,'fontsize',9.5);
        xtickangle(90);
        set(gcf,'color','w');
    end
end
    
end
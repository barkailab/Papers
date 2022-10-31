function [paraOrder, yTicks] = plotTrees(TFs, axesPosition, fig4,varargin)
% yTicks(1) = y position of the closer paralog!
load('/home/labs/barkailab/felixj/Documents/MATLAB/projects/Tamar/treeFiles.mat');
ip = inputParser;
ip.addParameter('colSpot', lines(2))
ip.parse(varargin{:})
colSpot = ip.Results.colSpot;

treeIdx = find(contains(treeFiles.name, TFs, 'IgnoreCase',true) & treeFiles.order == 3);
if numel(treeIdx) ==1
    currTree = treeFiles.tree{treeIdx};
    currTree=prune(currTree,find(~contains(get(currTree,'NodeNames'),{'Root','Branch','Klac','Tdel','Zrou','Scer','Suva','Skud','Smik','Node'})));
    plot(currTree)
    for ch =  get(gca,'Children')'
        if strcmp(class(ch),'matlab.graphics.primitive.Text')
            delete(ch)
        end
    end
    currTreePlot=copyobj(gca,fig4);
    close gcf
    set(currTreePlot,'Position',  axesPosition, 'YTick',[], 'XTick',[])
    axes(currTreePlot)
    hold on
    allNodes=get(currTree,'nodeName');
    allNodes=allNodes(~contains(allNodes,{'Root','Branch','Node'}));
    [~,typeIdx]=ismember(extractAfter(allNodes,'_'), extractAfter([treeFiles.Name1(treeIdx),treeFiles.Name2(treeIdx)],'_'));
    
    treeChildren = get(currTreePlot,'Children')';
    for chIdx=1:numel(treeChildren)
        ch = treeChildren(chIdx);
        if strcmp(class(ch),'matlab.graphics.chart.primitive.Line')
            if numel(ch.XData)==numel(allNodes) & numel(unique(diff(ch.YData)))==1
                for j=1:2
                    scatter(ch.XData(typeIdx==j),ch.YData(typeIdx==j),4,colSpot(j,:),'Marker',ch.Marker)
                    %text(mean(ch.XData(typeIdx==j)),mean(ch.YData(typeIdx==j)),protNames{j},'Color',colSpot(j,:),'FontSize',15,'HorizontalAlignment','center')
                    yTicks(j)=mean(ch.YData(typeIdx==j));
                    
                end
                scatter(ch.XData(startsWith(allNodes,'Klac')),ch.YData(startsWith(allNodes,'Klac')),4, 'ks')
            else
                ch.MarkerSize = ch.MarkerSize/4;
            end
            if ~strcmp(ch.MarkerFaceColor,'none')
                ch.MarkerFaceColor=.9*[1 1 1];
            end
            ch.Color=[.4 .4 .4];
            ch.MarkerEdgeColor=.5*[1 1 1];
        elseif strcmp(class(ch),'matlab.graphics.primitive.Text')
            delete(ch)
        end
    end
    if yTicks(1)>yTicks(2)
        set(currTreePlot, 'YDir','normal')
    end
    axis off; box off
    [~, paraOrder] =  ismember(extractAfter([treeFiles.Name1(treeIdx),treeFiles.Name2(treeIdx)],'_'), upper(TFs));
    xlim([-0.2 max(4.5, max(xlim))]);
else
    paraOrder = [1,2];
end
paraOrder(paraOrder==0) = setdiff([1,2], paraOrder);
end
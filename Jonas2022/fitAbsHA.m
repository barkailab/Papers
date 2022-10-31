function [yMatfit,yMatdT,a]=fitAbsHA(timeVec,yMat,varargin)
    ip=inputParser;
    ip.addParameter('figure',false)
    ip.parse(varargin{:});
    nClu=size(yMat,1)    ;
    lowB=[1 0.005 1 1.5 0];
    upperB=[1.1 10 inf 2 30];
    yMatfit=nan(size(yMat));
    yMatdT=nan(size(yMat));
    for i=1:nClu
        try
            [a{i},b(i)]=L5P(timeVec([1:end]),yMat(i,[1:end]),[],lowB,upperB);
            yMatfit(i,:)=a{i}(timeVec);
            yMatdT(i,:)=differentiate(a{i},timeVec);
        end
    end
    if ip.Results.figure
        figure
        lineCol=lines(4);
        c=0;
        for i=find(~cellfun('isempty',a))
            if c==12*4
                c=0;
                figure
            end
            c=c+1;
            subplot(3,4,ceil(c/4))
            plot(min(timeVec):max(timeVec),a{i}(min(timeVec):max(timeVec)),'Color',lineCol(mod(c-1,4)+1,:),'Linewidth',1)
            hold on
            scatter(timeVec,yMat(i,:),[],lineCol(mod(c-1,4)+1,:),'x')
            ylim([0.9 max([yMat(:);yMatfit(:)])])
            xlim(quantile(timeVec,[0 1]))
            xlabel('time in HU')
            ylabel(' abs DNA')
        end
    end
end
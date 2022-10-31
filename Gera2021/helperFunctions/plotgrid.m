function plotgrid(figure_data, varargin)
ip = inputParser;
ip.addParameter('color','k')
ip.parse(varargin{:})
hold on
    %along y axis
    for i = [0.5:1:size(figure_data,2)]
        plot( [i i],get(gca,'ylim'),'color', ip.Results.color,'linestyle','-');
    end
    %along x axis
    for i= [0.5:1:size(figure_data,1)]
        plot(get(gca,'xlim'), [i i],'color', ip.Results.color,'linestyle','-');
    end
end
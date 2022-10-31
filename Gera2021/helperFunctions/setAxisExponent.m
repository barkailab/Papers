function setAxisExponent(varargin)
ip = inputParser;
ip.addParameter('axes', gca);
ip.addParameter('axis', {'xy'});
ip.parse(varargin{:})
if contains(ip.Results.axis, 'x')
    magnitude = floor(log10(max(ip.Results.axes.XLim)));
    ip.Results.axes.XAxis.Exponent = magnitude;
end

if contains(ip.Results.axis, 'y')
    magnitude = floor(log10(max(ip.Results.axes.YLim)));
    ip.Results.axes.YAxis.Exponent = magnitude;
end

if contains(ip.Results.axis, 'c')
    magnitude = floor(log10(max(ip.Results.axes.Limits)));
    ip.Results.axes.Ruler.Exponent = magnitude;
end

end
function [xPoints,h0] = plotRaster_single(h0,MUFiring,fsamp,col,line_st,MUn)
ax = gca;    
spktimes = [];
spike_id = [];
for j=1:size(MUFiring,2)
    spktimes = [spktimes (MUFiring{1,j})'];
    spike_id = [spike_id j*ones(1,numel(MUFiring{1,j}))];
end

%plot
xPoints = [ spktimes; spktimes; NaN(size(spktimes)) ];
yPoints = [ MUn*spike_id-0.4;
            MUn*spike_id+0.4;
            NaN(size(spike_id)) ];
if ~isempty(fsamp)
    xPoints = (xPoints(:)-1)/fsamp;
else
    xPoints = xPoints(:);
end
yPoints = yPoints(:);
line(ax,xPoints,yPoints,'Color',col,'LineStyle',line_st);

end
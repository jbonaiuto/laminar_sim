function plot_fading_line(x, y, err, percentage, cm, line_style)

for idx=1:length(x)
    cidx=round(percentage(idx)*9)+1;
    color=cm(cidx,:);
    h=errorbar(x(idx),y(idx),err(idx),'color',color,'LineStyle','none','Marker','none','MarkerSize',6);
    if idx==1
        plot([x(idx) x(idx)+((x(idx+1)-x(idx))*0.5)], [y(idx) y(idx)+((y(idx+1)-y(idx))*0.5)], line_style, 'color', color);
    elseif idx==length(x)
        plot([x(idx-1)+((x(idx)-x(idx-1))*0.5) x(idx)], [y(idx-1)+((y(idx)-y(idx-1))*0.5) y(idx)], line_style, 'color', color);
    else
        plot([x(idx-1)+((x(idx)-x(idx-1))*0.5) x(idx)], [y(idx-1)+((y(idx)-y(idx-1))*0.5) y(idx)], line_style, 'color', color);
        plot([x(idx) x(idx)+((x(idx+1)-x(idx))*0.5)], [y(idx) y(idx)+((y(idx+1)-y(idx))*0.5)], line_style, 'color', color);
    end
end
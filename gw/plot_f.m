function [] = plot_f(name, xax, yax, tick_size)
%
%title(name);
xlabel(xax,'interpreter','latex');
ylabel(yax,'interpreter','latex');
set(gca,'TickLabelInterpreter','latex');
%ylh = get(gca,'ylabel');
%gyl = get(ylh);                                                         % Object Information
%ylp = get(ylh, 'Position');
%set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
set(gca,'FontSize',tick_size);
%set(gca, 'FontName', 'Times');
grid on;

end


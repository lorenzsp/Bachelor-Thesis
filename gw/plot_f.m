function [] = plot_f(name, xax, yax, tick_size)
%
%title(name);
xlabel(xax);
ylabel(yax);
ylh = get(gca,'ylabel');
%gyl = get(ylh);                                                         % Object Information
%ylp = get(ylh, 'Position');
%set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize',tick_size);
%set(gca, 'FontName', 'Times');
grid on;

end


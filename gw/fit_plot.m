function [fig_nums] = fit_plot(x,y,dy,fit_struct,varargin)
% function [fig_nums] = fit_plot(x,y,dy,fit_struct,varargin)
%
% funzione per produrre un grafico dei dati e del modello di un fit, 
% con i suoi residui
% input:
%   x,y: dati (indipendenti e dipendenti) oggetti del fit ad una retta
%   dy:  barre di incertezza nella variabile dipendente
%   fit_struct: structure output da "regressione_lineare"
% output: fig_nums vettore con [figure_num axes_dati axes_residui]
%
%   varargin (optional):   
%       'col'   prossimo argomento è il colore da usare
%       'dx'    barre di errore nella variabile indipendente
% 

col = 'b';
if length(varargin)>0
    for jj=1:length(varargin)
        if strcmp(varargin{jj},'dx')
            dx = varargin{jj+1};
        end
        if strcmp(varargin{jj},'col')
            col = varargin{jj+1};
        end
    end
end

if length(dy)==1
    dy = dy*ones(size(x));
end

fig=figure;
ax_top = axes('position',[.15 .45 .75 .45]);
if exist('dx')
    if length(dx) == 1
        dx = dx * ones(size(x));
    end
    errorbarxy(x,y,dx,dy, ... 
        'marker','.','linestyle','none','color',col);
else 
    errorbar(x,y,dy, ... 
        'marker','.','linestyle','none','color',col);
end
hold on;

% generate line from model
DX = max(x) - min(x);
x_mod = min(x)-0.05*DX + 1.1*DX*[0 1]';
y_mod = fit_struct.b + fit_struct.m*x_mod;
plot(x_mod,y_mod,'k');
ylabel('dati');
grid on;
xlim([min(x)-0.06*DX max(x)+0.06*DX]);
set(gca,'xticklabel','');


% generare residui ... 
res = y - (fit_struct.b + fit_struct.m*x);


% ricalcolare barre di incertezza per gli errori x (se presenti)
if exist('dx')
    dy = sqrt ( dy.^2 + (fit_struct.m*dx).^2);
end

% plot residui
ax_bot = axes('position',[.15 .15 .75 .23]);
errorbar(x,res,dy,...
        'marker','.','linestyle','none','color',col);
ylabel('residui');
xlabel('x');
grid on;
xlim([min(x)-0.06*DX max(x)+0.06*DX]);

fig_nums=[fig ax_top ax_bot];
%% wave polarizations animated
%angle parameter of the ring
theta = 0:0.03:2*pi;
%frequency of the wave
omega = 0.1;
% t time


figure();

n=1;
ii=1;
for ot=0:0.1:2*pi;
    % plus polarization
    h_plus = 0*cos(ot);
    h_times = 0.5*cos(ot);
    X = cos(theta) .* (1 + 0.5.*h_plus) + sin(theta).*(0.5.*h_times);
    Y = sin(theta) .* (1 - 0.5.*h_plus) + cos(theta).*(0.5.*h_times);

    CM = jet(150); % n+10 
    r=1;
    %for r=0.05:0.05:2
    plot(r.*X,r.*Y,'.','color',CM(n,:));%,ii.*ones(size(X))
    %hold on;
    %end
    grid on;
    %set(gca,'Visible','off')
    pbaspect([1 1 1]);
    %set(gca,'zticklabel',[])
    %set(gca,'xticklabel',[])
    %set(gca,'yticklabel',[])
    plot_f('','x','y',20);
    %xlabel('x');
    %ylabel('y');
    zlabel('t');
    %set(gca,'xtick',[])
    %set(gca,'ytick',[])
    %set(gca,'ztick',[])
    xlim([-2 2]);
    ylim([-2 2]);
    %pause(0.01)
    if ot==0
        gif('times_pol.gif','DelayTime',0.2,'LoopCount',5,'frame',gcf);
    else
    
    gif
    end

    %hold on;
    %n=n+1;
    clf;
    ii = ii + 0.2;
end





%% binary sources and GW
[C,ia,ib] = intersect(b4_p(:,1),b4_h(:,1));
t=C;
for n=1:10:length(C)
    
    subplot(2,1,1),
    plot(b4_p(ia(n),2),b4_p(ia(n),3),'.','MarkerSize',15);
    hold on;
    plot(-b4_p(ia(n),2),-b4_p(ia(n),3),'.','MarkerSize',15);
    axis equal
    xlim([-11 11]);
    ylim([-11 11]);
    
    subplot(2,1,2), 
    hold on;
    plot(b4_h(:,1),b4_h(:,2),'r');
    plot(b4_h(:,1),b4_h(:,3),'b');
    plot(t(n),b4_h(ib(n),2),'or','MarkerSize',4,'MarkerFaceColor','m') ;   
    plot(t(n),b4_h(ib(n),3),'or','MarkerSize',4,'MarkerFaceColor','m');
    xlim([0 400]);
    
    pause(0.01)
    clf;

end

for n=ib(end):24:24849
    
    subplot(2,1,1),
    plot(0,0,'.','MarkerSize',15);
    hold on;
    %plot(-b4_p(ia(n),2),-b4_p(ia(n),3),'.','MarkerSize',15);
    axis equal
    xlim([-11 11]);
    ylim([-11 11]);
    
    subplot(2,1,2), 
    hold on;
    plot(b4_h(:,1),b4_h(:,2),'r');
    plot(b4_h(:,1),b4_h(:,3),'b');
    plot(b4_h(n,1),b4_h(n,2),'or','MarkerSize',4,'MarkerFaceColor','m') ;   
    plot(b4_h(n,1),b4_h(n,3),'or','MarkerSize',4,'MarkerFaceColor','m');
    xlim([0 400]);
    
    pause(0.01)
    clf;

end

%%




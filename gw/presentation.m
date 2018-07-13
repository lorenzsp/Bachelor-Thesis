%% to create frames from gif
% gifsicle --unoptimize ../times_pol.gif | convert - frame-%d.png
%% wave polarizations animated
%angle parameter of the ring
theta = 0:0.1:2*pi;
%frequency of the wave
omega = 0.1;
% t time


figure();

n=1;
ii=1;
for ot=0:0.1:2*pi;
    % plus polarization
    h_plus = 0.5*cos(ot);
    h_times = 0.5*cos(ot-pi/2);
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
    pause(0.01)
    if ot==0
        gif('r_pol.gif','DelayTime',0.2,'LoopCount',5,'frame',gcf);
    else
    
    gif
    end

    %hold on;
    %n=n+1;
    clf;
    ii = ii + 0.2;
end





%% binary sources and GW
theta = 0:0.1:2*pi;
[C,ia,ib] = intersect(b3_p(:,1),b3_h(:,1));
t=C;
for n=1:10:length(C)
    
    subplot(2,2,1),
    plot(b3_p(ia(n),2),b3_p(ia(n),3),'.','MarkerSize',15);
    hold on;
    plot(-b3_p(ia(n),2),-b3_p(ia(n),3),'.','MarkerSize',15);
    %axis equal
        pbaspect([1 1 1]);

    xlim([-4 4]);
    ylim([-4 4]);
    plot_f('','$$x \; [M_{\odot}]$$','$$y \;[M_{\odot}]$$',16)
    grid on;
    
    subplot(2,2,2),
    % plus polarization
    h_plus = 200.*b3_h(ib(n),2);
    h_times = 200.*b3_h(ib(n),3);
    X = cos(theta) .* (1 + 0.5.*h_plus) + sin(theta).*(0.5.*h_times);
    Y = sin(theta) .* (1 - 0.5.*h_plus) + cos(theta).*(0.5.*h_times);
    plot(X,Y,'.');
    xlim([-2 2]);
    ylim([-2 2]);
        grid on;
            pbaspect([1 1 1]);


    
    
    subplot(2,2,[3,4]), 
    hold on;
    plot(b3_h(:,1),b3_h(:,2).*110/(100.*Mparsec),'r');
    plot(b3_h(:,1),b3_h(:,3).*110/(100.*Mparsec),'b');
    plot(t(n),b3_h(ib(n),2).*110/(100.*Mparsec),'or','MarkerSize',4,'MarkerFaceColor','m') ;   
    plot(t(n),b3_h(ib(n),3).*110/(100.*Mparsec),'or','MarkerSize',4,'MarkerFaceColor','m');
    ylim([-1.5 1.5].*1e-22)
    xlim([0 200]);
        grid on;
plot_f('','$$t \,[M_{\odot}]$$','$$ h(t,r=100 Mpc) $$',16)
    s= {['$$h_+$$'],...
    ['$$h_{\times}$$']};
    legend_f(s);
    
    %pause(0.01)
    if n==1
        gif('bbh-b3.gif','DelayTime',0.2,'LoopCount',5,'frame',gcf);
    else
    
    gif
    end
    clf;

end

for n=ib(end):40:length(b3_h(:,1))
    
    subplot(2,2,1),
    plot(0,0,'k.','MarkerSize',15);
    hold on;
    plot_f('','$$x \; [M_{\odot}]$$','$$y \;[M_{\odot}]$$',16)
xlim([-4 4]);
ylim([-4 4]);
    pbaspect([1 1 1]);

            grid on;

    
    subplot(2,2,2),
    h_plus = 200.*b3_h((n),2);
    h_times = 200.*b3_h((n),3);
    X = cos(theta) .* (1 + 0.5.*h_plus) + sin(theta).*(0.5.*h_times);
    Y = sin(theta) .* (1 - 0.5.*h_plus) + cos(theta).*(0.5.*h_times);
    plot(X,Y,'.');
    xlim([-2 2]);
    ylim([-2 2]);
        pbaspect([1 1 1]);
        grid on;
    
    subplot(2,2,[3,4]), 
    hold on;
    plot(b3_h(:,1),b3_h(:,2).*110/(100.*Mparsec),'r');
    plot(b3_h(:,1),b3_h(:,3).*110/(100.*Mparsec),'b');
    plot(b3_h(n,1),b3_h(n,2).*110/(100.*Mparsec),'or','MarkerSize',4,'MarkerFaceColor','m') ;   
    plot(b3_h(n,1),b3_h(n,3).*110/(100.*Mparsec),'or','MarkerSize',4,'MarkerFaceColor','m');
        ylim([-1.5 1.5].*1e-22)
    xlim([0 200]);
            grid on;
    plot_f('','$$t \,[M_{\odot}]$$','$$ h(t,r=100 Mpc) $$',16)
    s= {['$$h_+$$'],...
    ['$$h_{\times}$$']};
    legend_f(s);
    %115*/(100.*Mparsec)
    
    %pause(0.01)
    gif
    clf;

end
close all
%% travelling GW
%angle parameter of the ring
theta = 0:0.01:2*pi;
%frequency of the wave
omega = 3;
% t time

n=1;
ii=1;
z=0.5:0.07:15;
figure();
for t=0.01:0.05:2.0944
    % plus polarization
    

    plot3(0,2*cos(omega.*(t)),2*sin(omega.*(t)),'.','MarkerSize',20);
    hold on;
    plot3(0,-2*cos(omega.*(t)),-2*sin(omega.*(t)),'.','MarkerSize',20);
    axis equal
    CM = jet(250); % n+10 
    %for z=0.5:0.1:4;
    h_plus = -cos(2*omega.*(t-z ))./(z);
    h_times = -sin(2*omega.*(t-z ))./(z);
    %X = cos(theta) .* (1 + 0.5.*h_plus) + sin(theta).*(0.5.*h_times);
    %Y = sin(theta) .* (1 - 0.5.*h_plus) + cos(theta).*(0.5.*h_times);
    plot3(z,h_plus,0.*ones(size(z)));%,'.','color',CM(n,:));%z.*round(10*z)
    area(z,h_plus);

    plot3(z,0.*ones(size(z)),h_times); 
    fill3(z([1 1:end 1]),0.*ones(length(z)+2), [0 h_times 0], 'r');
    
    grid on;
    %set(gca,'Visible','off')
    pbaspect([1 1 1]);
    %set(gca,'zticklabel',[])
    %set(gca,'xticklabel',[])
    %set(gca,'yticklabel',[])
    %plot_f('','z','h',20);
    xlabel('z');
    ylabel('h_+');
    zlabel('h_\times');
    %set(gca,'xtick',[])
    %set(gca,'ytick',[])
    %set(gca,'ztick',[])
    xlim([0 6]);
    zlim([-2 2]);
    ylim([-2 2]);
    view(37.5000,30)
    %pause(0.00005)
    n=n+1;
    hold on;
    if t==0.01
        gif('gw_travelling.gif','DelayTime',0.2,'LoopCount',5,'frame',gcf);
    else
    
    gif
    end
    
    clf;
    ii = ii + 0.2;
end

%% BBH-b6 BBH-b7 BBH-b10


%% binary sources and radius
theta = 0:0.1:2*pi;
[C,ia,ib] = intersect(b10_p(:,1),b10_h(:,1));
t=C;
for n=1:50:length(C)
    
   
    
    subplot(2,1,1),
    % plus polarization
    plot(b10_p(:,1),b10_p(:,4));
    hold on
   plot(t(n),b10_p(n,4),'or','MarkerSize',4,'MarkerFaceColor','m') ;   
   plot_f('','$$t \; [M_{\odot}]$$','$$R$$',16)
    %xlim([0 1200]);
    %ylim([-2 2]);
        grid on;
     %       pbaspect([1 1 1]);


    
    
    subplot(2,1,2), 
    hold on;
    plot(b10_h(:,1),b10_h(:,2).*115/(100.*Mparsec),'r');
    plot(b10_h(:,1),b10_h(:,3).*115/(100.*Mparsec),'b');
    plot(t(n),b10_h(ib(n),2).*115/(100.*Mparsec),'or','MarkerSize',4,'MarkerFaceColor','m') ;   
    plot(t(n),b10_h(ib(n),3).*115/(100.*Mparsec),'or','MarkerSize',4,'MarkerFaceColor','m');
    ylim([-4 4].*1e-23)
    %xlim([0 1200]);
        grid on;
plot_f('','$$t \,[M_{\odot}]$$','$$ h(t,r=100 Mpc) $$',16)
    s= {['$$h_+$$'],...
    ['$$h_{\times}$$']};
    legend_f(s);
    
    pause(0.01)
    %if n==1
    %    gif('bbh-b3.gif','DelayTime',0.2,'LoopCount',5,'frame',gcf);
    %else
    %
    %gif
    %end
    clf;

end

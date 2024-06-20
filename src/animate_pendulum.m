function animate_pendulum(t,thet)

L = 1;
x = L*sin(thet);
y = -L*cos(thet);

figure
hold on
box on
axis square
set(gca,'XTick',[], 'YTick',[], ...
        'XLim',[-1,+1]*1.1, 'YLim',[-1,+1]*1.1)


plh_bar = plot([0,x(1)],[0,y(1)],'-b','LineWidth',1);   % Pendulum arm
plh_bob = plot(x(1),y(1),'.b','MarkerSize',25);         % Pendulum bob
plot(0,0,'.r','MarkerSize',15)                          % Pivot
drawnow

tic
while toc < t(end)
    [~,idx] = min(abs(t - toc));
    set(plh_bar,'XData',[0,x(idx)],'YData',[0,y(idx)])
    set(plh_bob,'XData',x(idx),'YData',y(idx))
    drawnow
end

set(plh_bar,'XData',[0,x(end)],'YData',[0,y(end)])
set(plh_bob,'XData',x(end),'YData',y(end))
drawnow
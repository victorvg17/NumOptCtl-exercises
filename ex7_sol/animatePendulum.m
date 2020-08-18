function [] = animatePendulum(Theta, dt, filename)
% Animate the pendulum trajectory, creates a gif
% input:
%     Theta:    trajectory of pendulum angles
%     dt:       time step between succesive values of theta (in animation)
%     filename: filename of created gif (string)


first_it = true;
for th = Theta
    figure(100); clf; hold on;

    plot( [0, -sin(th)], [0, cos(th)], 'ko-')
    plot(-sin(th), cos(th), 'o', 'MarkerSize', 15, 'MarkerFaceColor', 'b' )
    axis(1.5 * [-1, 1, -1, 1])
    set(gca,'XColor','none','YColor','none','TickDir','out')
    
    pause(dt)


    frame = getframe(gcf);

    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if first_it
      first_it = false;
      imwrite(imind,cm,filename,'gif','DelayTime', dt, 'Loopcount',inf);
    else
      imwrite(imind,cm,filename,'gif','DelayTime', dt, 'WriteMode','append');
    end
  

end


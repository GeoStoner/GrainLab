% Simple model of ripple creation in sand bend using constant flight
% distance.
% Ryan Stoner. March 30, 2016 for modeling in the Earth Sciences

clear
clc
%% Initialize
% Distances x, spacing dx, maximum and minimum distances xmin, xmax
% [L], meters

xmin = 0;
xmax = 500;
dx = 1;

% Setting up initial distance and height
x_axis = xmin:dx:xmax;                       
y_axis = ones(1,length(x_axis))+(rand(1,length(x_axis))-0.9); 

% Setting up impact trajectories as straight lines with slope m [L/L]. We
% also need to calculate the bmax, the maximum y intercept
m = - 0.2;
bmax = y_axis(length(y_axis))-m*x_axis(length(y_axis));

% Maximum (farthest impact) and minumum impact. For initial plotting
% purposes only
ymax = m*x_axis + bmax;
ymin = m*x_axis + y_axis(1);           

% Plot of initial conditions
figure(1)
plot(x_axis,y_axis,'k-')
ylim([0 5])
hold on
plot(x_axis,ymax,'b-')
ylim([0 10])
figure(1)
plot(x_axis,ymin,'g-')
ylim([0 5])
hold off

%% Looping
% n - number of grains
n = 10000;

% Outside of loop define width of box (grwspl) and average number of grains
% moved (avegrloss) to obtain height.
% Initial speed in vertical (vertv) and horizontal (horv) for grains,
% [L/T], m/s. Gravitational const. (g) [L/T^2], m/s^2. Ball. trajectories.
avegrloss = 8;
grwspl = 10;
vertv = 4;
horv = 10;
g = 9.8;

% Add nplot - Number of plots, tplot - grains between plot
nplot = 100;
tplot = n/nplot;
nframe = 0;
len = length(x_axis);


for i= 1:n
    % First calculate bmax - maximum intercept. In loop because it can
    % change, bhit - random intercept to determine line of random impact.
    % yhit - equation of that line
    bmax = y_axis(length(y_axis))-m*x_axis(length(y_axis));
    bhit = (bmax-y_axis(1))*rand(1) + y_axis(1);
    yhit = m*x_axis + bhit;
    
    % Find location of impact
    indbangrange = find(yhit-y_axis<y_axis);
    indbang = min(indbangrange);
    
    % Find number of grains lost, randomized (capture some of the
    % stochastic nature of system), rounded to make whole no. Divide by
    % width of grains in box/width of splash for a loss of certain height.
    loss = abs(avegrloss*rand(1));
    lossh = loss/grwspl; 
    
    % Loss of material in impact
    y_axis(indbang) = y_axis(indbang) - lossh;
    
    % Calculate grain trajectory and x index of landing location. Since we
    % have wraparaound boundary condition I calculate another trajectory,
    % which is shifted by "1 graph."
    ytraj = y_axis(indbang) + (vertv * (x_axis-indbang)/horv)-...
    ((1/2)*(g*(x_axis-indbang).^2)/(horv^2));

    vertv2 = vertv - (g*(x_axis(len-indbang+1))/horv);
    ytraj2 = ytraj(len)+(vertv2*(x_axis)/horv)-((1/2)*...
        (g*(x_axis).^2)/(horv^2));
    
    
    indlandrange1 = find(ytraj-y_axis>-0.0001);
    indlandrange2 = find(ytraj2-y_axis>-0.0001);
    indland1 = max(indlandrange1);
    indland2 = max(indlandrange2);
    indland = min([indland1 indland2]);
    % one of the two arrays should always be empty

    % Add grains splashed at a certain distance, wrap around boundary
    if isempty(indland)==1
        y_axis(indbang)=y_axis(indbang)+lossh;
    else
    y_axis(indland)= y_axis(indland)+lossh;
    end

    
    
    
    if(rem(nframe,tplot)==0)
    

    figure(2)
    plot(x_axis,y_axis,'k-')
    hold on
    plot(x_axis,yhit,'k-.')
    ylim([-1 100])
    fill(x_axis,y_axis,'g')
    ylim([-1 100])
    xlabel('distance (m)')
    ylabel('height (m)')
    title('Saltation of Grains to create ripples')
    pause(0.02)
    hold off
    end
% Converting time to a string to print string out
% trounded = round(n(i),0);
% strtime = num2str(trounded);

   nframe = nframe+1;
   y_axis(1)=y_axis(len);
end



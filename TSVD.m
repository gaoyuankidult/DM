% Choose relative noise level in simulated noisy data
noiselevel = 0.0001;

% Decide the measurement inverval 
t_start_pt = -0.2;
t_start_end = 1.2;
s_start_pt = -10;
s_start_end = 10;

% Construct target image at desired resolution
N_t      = 1024;
N_s      = 1024;

t_points = linspace(t_start_pt,t_start_end,N_t)';
s_points = linspace(s_start_pt,s_start_end,N_s)';

% Construct measurement. The vector s contains the linear coordinates
% of the Radon transform in pixel coordinates.

f_m = foo_fun(t_points);

% Conbstruct forward map

a = ones(N_s,N_t);
a(:,[1,end]) = 0.5;
a = a *t_points(end)/N_t;
A_TS =-s_points*t_points';
A_EXP = exp(A_TS);
A = a.*A_EXP;

% Plot to check the measurement of original graph
figure(1)
plot(t_points,f_m)

% Plot of forward construction(i.e. laplace transform)
figure(2)
m = A*f_m+ noiselevel * randn([length(s_points),1]);
plot(s_points,m)

% Plot of original graph's reconstruction
figure(3)
f_recon =pinv(A)*m;
plot(t_points,f_recon);

% Conpute the sigular values of A
[U,S,V] = svd(A);
svals = diag(S);

% Parameters for formatting the plotted figures. Modify if the linewidth, 
% font size or point size in the plots are not optimal in the figures.

lwidth = 1;
thickline = 2;
fsize = 12;
msize = 4;


% From here, we start the truncated SVD routine.
n = length(svals);
ivec = linspace(1,17,17);

%Iterate through the sigular values
for i = ivec
    
    % Plot current singular vector
    p1 = figure(4);
    clf
    
    % Reconstruct from noisy data using only large singular values.
    ns               = i;
    [row,col]        = size(S.');
    Dplus            = sparse(row,col);
    Dplus(1:ns,1:ns) = diag(1./svals(1:ns));
    rec              = V*Dplus*U.'*m;
    
    % Calculate relative error and plot reconstructed graph
    relerr = round(norm(rec-f_m)/norm(f_m)*100);
    plot(t_points,rec,'r','linewidth',lwidth)
    title([num2str(relerr), '% errors']);
    
    % Set plotting parameters
    p2 = figure(5);
    axis([0 1 -.5 2])
    set(gca,'xtick',[0 1/2 1],'fontsize',fsize)
    set(gca,'ytick',[0  1 2],'fontsize',fsize)
    set(gca,'xticklabel',{})
    set(gca,'yticklabel',{})
    set(gca,'PlotBoxAspectRatio' , [2 1 1])
    box off
    title(['Reconstruction (red) using ', num2str(i), ' singular vectors'])

    % Plot the used and unused singular values
    clf
    semilogy([1:n],svals,'k','linewidth',lwidth);
    hold on
    semilogy([1:i],svals(1:i),'r','linewidth',thickline);
    xlim([1 n])
    axis square;
    title('Used (red) and unused (black) singular values');
    saveas(p1,sprintf('rec_%d',i),'jpg');
    saveas(p2,sprintf('no_%d',i),'jpg');
    drawnow
end

% This procedure evaluates the objective form 
%
%      f(x) = 1 for 0 <= t <= 1 else f(x) = 0 
%
% Arguments:
% x         Evaluation point
%
% Returns:
% result_array 		value of the discrepancy function at point x
%

function result_array = foo_fun(t)
    result_array = arrayfun(@compare,t);
end 

function result = compare(t_ele)
    if 0 <= t_ele && t_ele <= 1
        result = 1;
    else
        result = 0;
    end 
end 


close all
clear

% import data from comsol
A0 = importdata('periodic_control_nov24b_control.txt');



k = 3e-2;

x_aux = [55.5 56 58:3:82];%  
x = 55*ones(2*length(x_aux),1);
x(2:2:end) = x_aux;


T0 = zeros(size(A0,1)-3,size(A0,2));

aux = size(T0,1)/length(x);

for l=1:size(T0,1)
    T0(l,:) = [str2double(A0{l+3,1}) str2double(A0{l+3,2})];
end

% create matrix with values of sigma
sigma = zeros(aux,length(x),length(k));
t = T0(1:aux,1);
n = 1;

for l = 1:length(x)
    sigma(:,l,1) = T0(1+aux*(n-1):aux*n,2);
    n = n + 1;
end



% plot tension propagation over time
for m = 1:2:length(x)
    leg = [];
    figure
    set(gcf, 'Position',  [50, 50, 800, 400])
    hold on 
    
    for l = 1:length(k)
        plot(t,sigma(:,m,l),'linewidth',1)
        plot(t,sigma(:,m+1,l),'linewidth',2)
        set(gca,'FontSize',14)
        xlabel('time [s]')
        ylabel('\sigma [pN/{\mu}m]')
        leg = [leg;['k = ' num2str(k(l),'%.3f') '[{\mu}m^2] active ']];
        leg = [leg;['k = ' num2str(k(l),'%.3f') '[{\mu}m^2] passive']];
    end
    
    legend(leg)
    title(['x = ' num2str(x(m+1),'%.2f') '[{\mu}]'])
end
hold off



col_plot =  [0.5,0.5,0.5];

% get the value of delta
type = {'control'};
aux = zeros(length(x),length(k));
fit_a = zeros(length(k),1);
fit_b = zeros(length(k),1);
fit_c = zeros(length(k),1);
fit_a_CI = zeros(length(k),1);
fit_b_CI = zeros(length(k),1);
fit_c_CI = zeros(length(k),1);
fit_rmse = zeros(length(k),1);
fit_b_SE = zeros(length(k),1);
fit_a_SE = zeros(length(k),1);
fit_c_SE = zeros(length(k),1);
fit_df = zeros(length(k),1);
est = 0.95;
peaks = zeros(length(x)/2,length(k));
passive = zeros(length(x)/2,length(k));
p_x = x(2:2:end)-x(1);
x_plot = p_x(1):0.1:p_x(end);

exp_decay = 'a*exp(-x/b)+c';
startPoints = [1 1 0];


aux(:,:) = max(sigma(:,1:length(x),:));
figure
set(gcf, 'Position',  [50, 50, 400, 400])
hold on 
for j = 1:length(k)
    peaks(:,j) = (aux(1:2:end,j)-sigma(1,1:2:length(x),j)');
    passive(:,j) =(aux(2:2:end,j)-sigma(1,1:2:length(x),j)');
    p_y = 100*passive(:,j)./peaks(:,j);
    [f1,gof1] = fit(p_x,p_y,exp_decay,'Start',startPoints,'Weights',passive(:,j));
    fit_a(j,1) = f1.a;
    fit_b(j,1) = f1.b;
    fit_c(j,1) = f1.c;
    fit_rmse(j,1) = gof1.rmse;

    bb = confint(f1,est);
    fit_a_CI(j,1) = f1.a-bb(1,1);
    fit_b_CI(j,1) = f1.b-bb(1,2);
    fit_c_CI(j,1) = f1.c-bb(1,3);
   
    
    fit_df(j,1) = gof1.dfe;
    t = tinv((1+est)/2, fit_df(j,1)); 
    fit_b_SE(j,1) = (bb(2,2)-bb(1,2)) ./ (2*t);
    fit_a_SE(j,1) = (bb(2,1)-bb(1,1)) ./ (2*t);
    fit_c_SE(j,1) = (bb(2,3)-bb(1,3)) ./ (2*t);

    aux_size = (passive(:,j)<=3);
    if sum(aux_size)>0
        plot(p_x(aux_size),p_y(aux_size),'o','markerfacecolor',...
            col_plot(j,:),'MarkerSize',6,...
            'markeredgecolor',col_plot(j,:),...
            'DisplayName',strcat(type(j),", $k$ =" + k(j), ...
            "($\mu$m$^2$), $\delta$ =" + round(f1.b,1), "($\mu$m)"));
    end
    aux_size = (passive(:,j)>3 & passive(:,j)<=6);
    if sum(aux_size)>0
        h = plot(p_x(aux_size),p_y(aux_size),'o','markerfacecolor',...
            col_plot(j,:),'MarkerSize',9,...
            'markeredgecolor',col_plot(j,:));
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
    aux_size = (passive(:,j)>6 & passive(:,j)<=9);
    if sum(aux_size)>0
        h1 = plot(p_x(aux_size),p_y(aux_size),'o','markerfacecolor',...
            col_plot(j,:),'MarkerSize',12,...
            'markeredgecolor',col_plot(j,:));
        h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
    

    h2 = plot(x_plot,f1.a*exp(-x_plot./f1.b)+f1.c,'color',col_plot(j,:),...
        'LineWidth',1);
    h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
end

set(gca,'FontSize',16)
xlabel('Tether separation ({\mu}m)')
ylabel('Tension propagation (%)')
legend('Interpreter',"latex"); 
axis([0 27 0 75])








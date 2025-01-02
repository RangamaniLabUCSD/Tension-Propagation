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
n = 1;
t = T0(1:aux,1);

for l = 1:length(x)
    sigma(:,l,end) = T0(1+aux*(n-1):aux*n,2);
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


col_plot =  [0.5 0.5 0.5];
% get the value of alpha
type = {'control'};
aux = zeros(length(x),length(k));
fit_a = zeros(length(k),1);
fit_b = zeros(length(k),1);
fit_adjrsquare = zeros(length(k),1);
peaks = zeros(length(x)/2,length(k));
passive = zeros(length(x)/2,length(k));
p_x = x(2:2:end)-x(1);
x_plot = p_x(1):0.1:p_x(end);


aux(:,:) = max(sigma(:,1:length(x),:));
figure
set(gcf, 'Position',  [50, 50, 400, 400])
hold on 
for j = 1:length(k)
    peaks(:,j) = (aux(1:2:end,j)-sigma(1,1:2:length(x),j)');
    passive(:,j) =(aux(2:2:end,j)-sigma(1,1:2:length(x),j)');
    p_y = 100*passive(:,j)./peaks(:,j);
    [f1,gof1] = fit(p_x,p_y,'exp1','Weights',passive(:,j));
    fit_a(j,1) = f1.a;
    fit_b(j,1) = f1.b;
    fit_adjrsquare(j,1) = gof1.adjrsquare;

    aux_size = (passive(:,j)<=3);
    if sum(aux_size)>0
        plot(p_x(aux_size),p_y(aux_size),'o','markerfacecolor',...
            col_plot(j,:),'MarkerSize',6,...
            'markeredgecolor',col_plot(j,:),...
            'DisplayName',strcat(type(j),", ", "$\alpha$ =" + -1/f1.b, "($\mu$m)"));
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

    
    h2 = plot(x_plot,f1.a*exp(f1.b.*x_plot),'color',col_plot(j,:),...
        'LineWidth',1);
    h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
end

set(gca,'FontSize',16)
xlabel('dist ({\mu}m)')
ylabel('Peak \sigma_{passive}/Peak \sigma_{active} (%)')
legend('Interpreter',"latex"); 
axis([0 27 0 75])

mean(peaks)
std(peaks)







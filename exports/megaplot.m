Volume = 3.209553837;

%% load data
file1 = './team/fine_probes_p.csv';
file2 = './team/medium_probes_p.csv';
file3 = './team/coarse_probes_p.csv';

nCells1 = 712233;
nCells2 = 211032;
nCells3 = 62528;

tcol = 1;

data1 = readmatrix(file1);
data2 = readmatrix(file2);
data3 = readmatrix(file3);

h1 = (Volume/nCells1)^(1/3)
h2 = (Volume/nCells2)^(1/3)
h3 = (Volume/nCells3)^(1/3)

% always h_coarse / h_fine
r21 = h2/h1
r32 = h3/h2

tq = [1.05344, 1.77, 2.90667, 3.964, 4.864, 5.856, 6.85333];

p_mean_vec = zeros(8,1);
GCI_mat = zeros(length(tq),8);
Phi_mat = zeros(length(tq),8);

for col = 2:9

    % interpol
    
    t1 = data1(:,tcol);
    Phi1 = data1(:,col);
    
    t2 = data2(:,tcol);
    Phi2 = data2(:,col);
    
    t3 = data3(:,tcol);
    Phi3 = data3(:,col);
    
    Phi1 = interp1(t1, Phi1, tq, 'linear');
    Phi2 = interp1(t2, Phi2, tq, 'linear');
    Phi3 = interp1(t3, Phi3, tq, 'linear');

    Phi_mat(:,col-1) = Phi1;
    
    N = length(Phi1);
    p = zeros(N,1);
    e21_a = zeros(N,1);
    for i = 1:N
    
        a=(Phi3(i)-Phi2(i))/(Phi2(i)-Phi1(i));
        s=sign(a);
        f=@(x)x-abs(log(abs(a))+log((r21^x-s)/(r32^x-s)))/log(r21);
        p(i)=fzero(f,1);
        phi_ext_21=(r21^p(i)*Phi1(i)-Phi2(i))/(r21^p(i)-1);
        phi_ext_32=(r32^p(i)*Phi2(i)-Phi3(i))/(r32^p(i)-1);
        e21_a(i)=abs((Phi1(i)-Phi2(i))/Phi1(i));
        e_ext=abs((phi_ext_21-Phi1(i))/phi_ext_21);
    
    end

    p_mean = mean(p)
    p_mean_vec(col-1) = p_mean;
    GCI=1.25*e21_a/(r21.^p_mean-1);
    GCI_mat(:,col-1) = GCI;

end

p_mean_vec


%% meag plot
val_p = readtable("exp_data.txt");

fig = figure
tiledlayout(2,2)
% 2-5, 6-9
for col = 6:9
    nexttile
    col_name_val = val_p.Properties.VariableNames{col};

    plot(val_p.Time_s_, val_p.(col_name_val), '--', 'LineWidth', lwidth, 'Color', smag)
    hold on
    plot(data1(:,1), data1(:,col), 'LineWidth', lwidth, 'Color', sblau)
    errorbar(tq, Phi_mat(:,col-1), GCI_mat(:,col-1).*Phi_mat(:,col-1), 'Color', sblau, "LineStyle","none", 'LineWidth', 1)
    hold off

    title(sprintf('Probe %d',col-1), 'Interpreter', 'latex')
    xlabel('time / s', 'Interpreter', 'latex')
    ylabel('$p$ / Pa', 'Interpreter', 'latex')
    grid on
    set(gca, 'FontSize', 15)
end
lgd = legend('$p_\mathrm{Simulation}$', '$p_\mathrm{Experiment}$', 'Interpreter', 'latex');
lgd.Layout.Tile= "south";

%exportgraphics(fig, './figures/p_probes_1_4.pdf', ContentType='vector')
exportgraphics(fig, './figures/p_probes_5_8.pdf', ContentType='vector')
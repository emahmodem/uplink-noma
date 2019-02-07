function visualize_the_network(params)
simulation_area =  (params.simulation_area_side(2) - params.simulation_area_side(1))^2;

mu_b = params.LA_B * simulation_area;
mu_h = params.LA_H * simulation_area;
mu_m = params.LA_M * simulation_area;
%mu_c = params.LA_C * simulation_area;
N_cells = poissrnd(mu_b);
while(N_cells < 1)
    N_cells = poissrnd(mu_b);
end
N_users_HTC = poissrnd(mu_h);
N_users_MTC = round(params.rho_m * poissrnd(mu_m));
%N_nodes_CH = round(mu_c);    % Fixed Privilged Cluster Heads

locations.BS = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_cells, 2);
locations.HTC = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users_HTC, 2);
locations.MTC = params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_users_MTC, 2);
%locations.CH= params.simulation_area_side(1) + 2 * params.simulation_area_side(2).* rand(N_nodes_CH, 2);

[association,loads] = associate_nodes(locations,params.aggregation_mode);

%plot (0,0,'X','MarkerSize',20,'LineWidth',4,'Color','black')
hold on;

MTC = plot_nodes(locations.MTC(:,1),locations.MTC(:,2),'.m',30);
HTC = plot_nodes(locations.HTC(:,1),locations.HTC(:,2),'sg',15);
cells = plot_nodes(locations.BS(:,1),locations.BS(:,2),'^b',15);
active_cells = plot_nodes(locations.BS(association.Active_Servers,1),locations.BS(association.Active_Servers,2),'^r',15);
annotate_node_id(locations.BS,15)
annotate_node_id(locations.HTC,15)
annotate_node_id(locations.MTC,10)
if(strcmp(params.aggregation_mode,'AGGREGATION') )
    CH = plot_nodes(locations.CH(:,1),locations.CH(:,2),'*m',15);
    annotate_node_id(locations.CH,15)
    legend ({'MTC'  , 'HTC' ,'BS','CH' })
else
    legend ({'MTC'  , 'HTC' ,'BS','Active BS'})
end
connect_nodes_to_servers(locations,association,params.aggregation_mode)

axis equal;
X_lim = params.simulation_area_side;
Y_lim = params.simulation_area_side;
axis([X_lim Y_lim]);
xlabel('x (meters)');
ylabel('y (meters)');
set(gca, 'FontName', 'Arial');
set(gca, 'FontSize', 20);
set(gca, 'FontWeight', 'Bold');

grid on;

%figure;
%bar(loads.BS)

% xlabel('Cell ID');
% ylabel('Number of connected nodes');
% set(gca, 'FontName', 'Arial');
% set(gca, 'FontSize', 20);
% set(gca, 'FontWeight', 'Bold');
% legend ({'MTC nodes'  , 'HTC users' })
% grid on;
end
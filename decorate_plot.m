function decorate_plot(h,deco)
    dark_green = [0 0.5 0];
    set(h([3 4]),'Color' , dark_green)
    %set(h([2]),'Color' , dark_green)
    set(h,'MarkerSize',15);
    set(h,'LineWidth',4);
    set(gca, 'FontSize', 30);
    set(gca, 'FontWeight', 'Bold');
    set(gca, 'LineWidth', 2);
    set(gca, 'GridAlpha', 0.5);
    set(gca, 'MinorGridAlpha', 0.5);
    lgd = legend(deco.legend.labels,'FontSize',25,'FontWeight','bold','Location',deco.legend.location,'Interpreter','LaTex');
    lgd.FontSize = 18;
    xlabel(deco.xlabel,'Interpreter','LaTex');
    ylabel(deco.ylabel,'Interpreter','LaTex');
    title(deco.title,'Interpreter','LaTex')
    grid on;
end



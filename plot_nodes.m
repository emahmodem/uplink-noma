function y = plot_nodes(x,y,marker,size)
    v = plot(x,y,marker);
    set(v,'MarkerSize',size);
    set(v,'LineWidth',4);
end
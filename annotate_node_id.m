function annotate_node_id(nodes,fontsize)
    for i = 1:size(nodes,1)
        text(nodes(i,1) + 1,nodes(i,2) + 1,num2str(i),'Color','black','FontSize',fontsize);
    end
end

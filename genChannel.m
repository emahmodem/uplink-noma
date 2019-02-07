function H = genChannel(numCells,numUsers)
    H = exprnd(1,numCells,numUsers);
    save('H.txt','H','-ascii') ;
end


function  plotCdGraphs()
%Kyle Simonis, 
%Created: ~2/21/2013
%Last Updated: 4/22/13
%Repeatedly executes the function calculateCd to calculate the coefficient
%of drag v. M# for different values of sigmaN, sigmaT, etc., and plotting
%separate graphs for each.

willLoop = true;

while willLoop == true
    [diskPerp, diskParallel, sphere, sharpCone, circCyl, M] = calculateCd(false);

    %Plot the graph of cD v. M (with the exception of the disk with parallel flow, which
    %is in its own subplot
    subplot(2, 1, 1);
    plot(M, diskPerp, M, sphere, M, sharpCone, M, circCyl);
    
    subplot(2,1,2);
    plot(M, diskParallel);

    %Asks if there is different data to input. If so, loop back. Otherwise,
    %end.
    if (input('To Exit, Type 0: ') == 0)
        willLoop = false;
    end



end
end


module EconTools

using Lib



function ces(quantities, sigma, integrate = True)
#========================================================================================
ces - Calculates the utility of a CES utility function

Inputs:
    quantities - a N x 1 vector which corresponds to the quantities consumed



=########################################################################################

end


function frechet(z, T, theta)
    return exp(-T*z^(-theta))



end
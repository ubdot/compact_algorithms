function xd = denorm(x,lowLim,upLim)
%Denormalize the sampled value, the value returned is in the
%domain of the search problem
%   Detailed explanation goes here
    len = upLim-lowLim;
    xd = .5*x.*(len) + lowLim + len/2;
end
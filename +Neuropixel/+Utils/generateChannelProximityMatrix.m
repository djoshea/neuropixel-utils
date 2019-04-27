function [proxMat, distMat] = generateChannelProximityMatrix(channelMap, distanceThresh, channelIds)
% generate an N x N matrix which is true where two channels are within
% distanceThresh of each other (in microns). By default uses connected
% channels but other channel inds may be specified

if nargin < 3
    channelIds = channelMap.connectedChannels;
end

channelInds = channelMap.lookup_channelIds(channelIds);
x = channelMap.xcoords(channelInds);
y = channelMap.ycoords(channelInds);
N = numel(x);

X = repmat(x, 1, N);
Y = repmat(y, 1, N);

distSq = (X - X').^2 + (Y - Y').^2;
proxMat = distSq < distanceThresh.^2;

if nargout >= 2
    distMat = sqrt(distSq);
end

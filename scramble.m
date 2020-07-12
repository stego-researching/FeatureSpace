function [ stream, sed ] = scramble(stream, sed)
%   scramble is the function that scramble the pixel order
%   @para        stream: bit stream
%                      seed : is a random seed

n = size(stream,1);                    % use n pixels, n is the size of cover

if ~exist('sed','var')
    sed  = sum(100*clock);
end

RandStream.setGlobalStream(RandStream('mt19937ar','seed',sed));  % scramble  seed

p = randperm(n);                          % scramble  seed

stream = stream(p,:);                         % scramble

end


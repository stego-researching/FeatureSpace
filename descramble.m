function [ stream ] = descramble(stream, sed)
%   scramble is the function that scramble the pixel order
%   @para        stream: bit stream
%                      seed : is a random seed

n = size(stream,1);                    % use n pixels, n is the size of cover

RandStream.setGlobalStream(RandStream('mt19937ar','seed',sed));  % scramble  seed

p = randperm(n);                          % scramble  seed

[~,img_indx]  = sort(p);               % Descramble used same seed
    
stream = stream(img_indx,:);

end


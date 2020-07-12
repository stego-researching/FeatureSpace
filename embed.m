function [ stego ] = embed(cover, m, flag)
%EMBED Summary of this function goes here
%   Detailed explanation goes here

if exist('flag', 'var')
    if strcmp(lower(flag), 'multi')
        stego = multi_embed(cover, m);
        return;
    end
end

% load necessary parameters
load('para.mat');

% generate distortiom map according to the image
dmap = gen_distmap_single(single(cover), kern, MU_R, coeff, w, max_v, pattern_select, pattern_map);
dmap = (dmap * 13).^ 0.5 + 0.5;

% random messages
message  = uint8(rand(m, 1));

stego = single_embed(cover, message, dmap);

end


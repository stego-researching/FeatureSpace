function [ stego ] = multi_embed(cover, m)
%MULTI_EMBED Summary of this function goes here
%   Detailed explanation goes here

% load necessary parameters
load('para.mat');

stego = cover;
message  = uint8(rand(m/16, 16));

for X_start = 1:4
    for Y_start = 1:4
        
        % generate dmap: costs of flipping pixels
        dmap_tmp = gen_distmap_multi(single(stego), kern, MU_R, coeff, w, max_v, pattern_select,...
            pattern_map, single(X_start - 1), single(Y_start - 1));
        
        dmap_tmp = dmap_tmp(X_start:4:end, Y_start:4:end);
        dmap_tmp = (dmap_tmp * 13) .^ 0.5 + 0.5;
        
        % embed message using STC
        cover_part = stego(X_start:4:end, Y_start:4:end);
        stego_part = single_embed(cover_part, message(:, (X_start -1)*4 + Y_start), dmap_tmp);
        stego(X_start:4:end, Y_start:4:end) = stego_part;
    end
end

stego = logical(stego);

end


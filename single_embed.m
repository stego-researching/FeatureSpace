function [ stego ] = single_embed(cover, message, dmap)
%SINGLE_EMBED Summary of this function goes here
%   Detailed explanation goes here

coverStream = uint8(cover(:));
n = size(coverStream,1);
alpha = size(message,1)/n;
dmap = double(dmap(:));
h = 10;

% scramble
[coverStream, sed] = scramble(coverStream);
dmap = scramble(dmap, sed);

tic;
[dist, stego] = stc_embed(coverStream, message, dmap, h); % embed message

fprintf('distortion per cover element = %f\n', dist / n);
fprintf('        embedding efficiency = %f\n', alpha / (dist / n));
fprintf('                  throughput = %1.1f Kbits/sec\n', n / toc() / 1024);

message2 = stc_extract(stego, size(message,1), h); % extract message
if all(message == message2)
    disp('Message has been extracted correctly.');
else
    error('Some error occured in the extraction process.');
end

% descramble
stego  = descramble(stego, sed);

stego = logical(reshape(stego, size(cover)));

end


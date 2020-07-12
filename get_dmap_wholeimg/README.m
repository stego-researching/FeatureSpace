clear;clc;
load('input_para.mat');
img = imread('1.bmp');
dmap = gen_distmapC_FS(single(img), single(kern), single(mu), single(coeff)...
    , single(w), single(max_v), single(pattern_select), single(pattern_map));
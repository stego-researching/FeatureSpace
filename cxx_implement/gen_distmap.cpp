#include <assert.h>
#include <memory.h>
#include <time.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
#include "mex.h"

using namespace std;

typedef unsigned int u32;
typedef unsigned short u16;
typedef unsigned char u8;

#define MATLAB_IN_IMG 0
#define MATLAB_IN_KERN 1
#define MATLAB_IN_MU 2
#define MATLAB_IN_COEFF 3
#define MATLAB_IN_W 4
#define MATLAB_IN_MAX_V 5
#define MATLAB_IN_SELECT 6
#define MATLAB_IN_MAPPER 7

#define MATLAB_OUT_DMAP 0

void print_vector(vector<float> vec, u32 start, u32 end) {
    for (u32 i = start; i < end; i++) {
        cout << vec[i] << ", ";
    }
    cout << endl;
}

void print_vector(vector<vector<float>> vec, u32 startx, u32 starty, u32 endx,
                  u32 endy) {
    for (u32 i = startx; i < endx; i++) {
        for (u32 j = starty; j < endy; j++) {
            cout << vec[i][j] << ", ";
        }
        cout << endl;
    }
    cout << endl;
}

/** Generate the "padded" image of size
(height+vector<vector<float>>*2)*(width+vector<vector<float>>*2).
@param im vector<vector<float>> within which image data is stored */
void padding(vector<vector<float>>& img, u32 edge_height, u32 edge_width) {
    vector<vector<float>> tmp(img);
    img.clear();
    img.resize(tmp.size() + 2 * edge_height);
    u32 resize_inside = tmp[0].size() + 2 * edge_width;
    for (u32 k = 0; k < img.size(); ++k) {
        img[k].resize(resize_inside);
    }
    /*now, allocate the space for a local "padded" copy of the image and copy
    it there, and padding value is '0' */
    for (u32 i = 0; i < tmp.size(); i++) {
        for (u32 j = 0; j < tmp[i].size(); j++) {
            img[i + edge_height][j + edge_width] = tmp[i][j];
        }
    }
}

// get the norm of the vector
vector<float> get_norm_value(vector<vector<float>>& w) {
    vector<float> norm_value(w.size(), 0);
    for (u32 i = 0; i < w.size(); i++) {
        float norm_w = 0;
        for (u32 j = 0; j < w[i].size(); j++) {
            norm_w = norm_w + (w[i][j] * w[i][j]);
        }
        norm_value[i] = sqrt(norm_w);
    }

    return norm_value;
}

float transfer2frcy(vector<float>& vec, vector<float>& max_v) {
    float sum_row = 0;
    for (u32 i = 0; i < vec.size(); i++) {
        sum_row += vec[i];
    }

    for (u32 i = 0; i < vec.size(); i++) {
        vec[i] = (vec[i] / sum_row) / max_v[i];
    }

    return sum_row;
}

u32 get_kern_offset(vector<vector<float>>& kern, bool flag) {
    u32 offset = 0;
    if ((flag ? kern[0].size() : kern.size()) % 2 == 0) {
        offset = 1;
    }
    return offset;
}

// get convlution result, img and kern
vector<vector<float>> imfilter(vector<vector<float>> img,
                               vector<vector<float>>& kern) {
    vector<vector<float>> indx(img.size(), vector<float>(img[0].size(), 0));

    if (img.size() < kern.size() + 2 || img[0].size() < kern[0].size() + 2) {
        cout << "Error: img should have bigger size, our imfilter can not "
                "handle the small img."
             << endl;
        exit(0);
    }

    u32 edge_height = (u32)(kern.size() / 2);
    u32 edge_width = (u32)(kern[0].size() / 2);

    // use padding function
    padding(img, edge_height, edge_width);

    u32 offset_height = get_kern_offset(kern, false);
    u32 offset_width = get_kern_offset(kern, true);

    // get the convlution result
    for (u32 m = 0; m < indx.size(); m++) {
        for (u32 n = 0; n < indx[m].size(); n++) {
            float conv_indx = 0;
            for (u32 i = 0; i < kern.size(); i++) {
                for (u32 j = 0; j < kern[i].size(); j++) {
                    conv_indx +=
                        kern[i][j] *
                        img[m + i + offset_height][n + j + offset_width];
                }
            }
            indx[m][n] = conv_indx;
        }
    }

    return indx;
}

vector<float> get_pattern_num(vector<vector<float>> img,
                              vector<vector<float>>& kern) {
    u32 count = 0;
    for (u32 i = 0; i < kern.size(); i++) {
        for (u32 j = 0; j < kern[i].size(); j++) {
            if (kern[i][j] != 0) {
                count++;
            }
        }
    }

    vector<float> pattern_num((u32)pow(2, count));

    vector<vector<float>> img_indx = imfilter(img, kern);

    u32 img_hei_bound = (u32)(kern.size() / 2);
    u32 img_wid_bound = (u32)(kern[0].size() / 2);

    for (u32 m = img_hei_bound; m < img.size() - img_hei_bound; m++) {
        for (u32 n = img_wid_bound; n < img[m].size() - img_wid_bound; n++) {
            u32 pos = (u32)img_indx[m][n];
            pattern_num[pos] = pattern_num[pos] + 1;
        }
    }

    return pattern_num;
}

// remove the idex of pattern which is not need
void PCA_mapper(vector<float>& indx_num, vector<float>& mu,
                vector<vector<float>>& coeff) {
    for (u32 i = 0; i < indx_num.size(); i++) {
        indx_num[i] = indx_num[i] - mu[i];
    }

    vector<float> indx_num_cp(indx_num);
    indx_num.clear();
    indx_num.resize(coeff[0].size());
    for (u32 j = 0; j < coeff[0].size(); j++) {
        float tmp = 0;
        for (u32 i = 0; i < coeff.size(); i++) {
            tmp += indx_num_cp[i] * coeff[i][j];
        }
        indx_num[j] = tmp;
    }

    indx_num_cp.clear();
}

// remove the idex of pattern which is not need
void PCA_mapper(vector<float>& indx_num_query, vector<float>& indx_num,
                vector<float>& indx_num_af, vector<float>& change_indx,
                vector<vector<float>>& coeff, vector<float>& max_v,
                float row_sum) {
    vector<float> indx_num_af_cp(indx_num_af);
    indx_num_af.clear();
    indx_num_af.resize(coeff[0].size());
    for (u32 j = 0; j < indx_num_af.size(); j++) {
        indx_num_af[j] = indx_num_query[j];
        for (u32 i = 0; i < change_indx.size(); i++) {
            indx_num_af[j] =
                indx_num_af[j] +
                (indx_num_af_cp[change_indx[i]] - indx_num[change_indx[i]]) /
                    row_sum / max_v[change_indx[i]] * coeff[change_indx[i]][j];
        }
    }
}

vector<float> mat2vec(vector<vector<float>>& matrix) {
    vector<float> vec(matrix.size() * matrix[0].size());
    u32 k = 0;
    for (u32 i = 0; i < matrix.size(); i++) {
        for (u32 j = 0; j < matrix[i].size(); j++) {
            vec[k++] = matrix[i][j];
        }
    }
    return vec;
}

vector<float> get_distance(vector<float>& indx_num, vector<vector<float>>& w,
                           vector<float>& norm_vector) {
    vector<float> distance(w.size());
    for (u32 i = 0; i < distance.size(); i++) {
        for (u32 j = 0; j < indx_num.size(); j++) {
            distance[i] += indx_num[j] * w[i][j];
        }

        distance[i] = distance[i] / norm_vector[i];
    }

    return distance;
}

void pattern_select(vector<float>& indx, vector<float>& select) {
    vector<float> indx_cp(indx);
    indx.clear();
    indx.resize(select.size());

    for (u32 i = 0; i < indx.size(); i++) {
        indx[i] = indx_cp[(u32)select[i]];
    }
    indx_cp.clear();
}

void mapper(vector<float>& indx, vector<float>& pmap) {
    for (u32 i = 0; i < indx.size(); i++) {
        indx[i] = pmap[(u32)indx[i]];
    }
}

void get_change_indx(vector<float>& indx_be, vector<float>& indx_af,
                     vector<float>& result) {
    for (u32 i = 0; i < indx_be.size(); i++) {
        result.push_back(indx_be[i]);
    }

    for (u32 i = 0; i < indx_af.size(); i++) {
        result.push_back(indx_af[i]);
    }

    sort(result.begin(), result.end());

    for (u32 i = 0; i < result.size() - 1; i++) {
        if (result[i] == result[i + 1]) {
            result.erase(result.begin() + i + 1);
            i--;
        }
    }

    if (result[0] == -1) {
        result.erase(result.begin());
    }
}

void gen_distortion_map(vector<vector<float>>& dmap, vector<vector<float>>& img,
                        vector<vector<float>>& kern, vector<float>& mu,
                        vector<vector<float>>& coeff, vector<vector<float>>& w,
                        vector<float>& max_v, vector<float>& select,
                        vector<float>& pmap) {
    u32 edge_height = (u32)(kern.size() / 2);
    u32 edge_width = (u32)(kern[0].size() / 2);

    padding(img, edge_height, edge_width);

    // Get the indx_num and calculate the distance_be before flipping
    vector<float> norm_vector = get_norm_value(w);

    vector<float> indx_num = get_pattern_num(img, kern);
    pattern_select(indx_num, select);

    vector<float> indx_num_be(indx_num);
    float row_sum = transfer2frcy(indx_num_be, max_v);
    PCA_mapper(indx_num_be, mu, coeff);

    vector<float> indx_num_query(indx_num_be);

    // get distance
    vector<float> distance_be = get_distance(indx_num_be, w, norm_vector);

    vector<vector<float>> indx = imfilter(img, kern);

    // set kern to be a vector
    vector<float> kern_vec = mat2vec(kern);

    u32 offset_height = get_kern_offset(kern, false);
    u32 offset_width = get_kern_offset(kern, true);

    vector<float> indx_be(kern.size() * kern[0].size());
    vector<float> indx_af(kern.size() * kern[0].size());

    // start generating the distortion scores of every pixel
    for (u32 i = edge_height; i < img.size() - edge_height; i++) {
        for (u32 j = edge_width; j < img[i].size() - edge_width; j++) {
            // get the indx_be from indx, coressponding to i, j
            u32 tmp_x = 0, tmp_y = 0;
            for (u32 k = 0; k < indx_be.size(); k++) {
                indx_be[k] = indx[i + edge_height - offset_height - tmp_x]
                                 [j + edge_width - offset_width - tmp_y++];
                if ((k + 1) % kern[0].size() == 0) {
                    tmp_y = 0;
                    tmp_x++;
                }
            }

            // get indx_af from the relationship bewteen indx_be and kern
            if (img[i][j] == 1) {
                for (u32 k = 0; k < indx_af.size(); k++) {
                    indx_af[k] = indx_be[k] - kern_vec[k];
                }
            } else {
                for (u32 k = 0; k < indx_af.size(); k++) {
                    indx_af[k] = indx_be[k] + kern_vec[k];
                }
            }

            mapper(indx_be, pmap);
            mapper(indx_af, pmap);

            // Change the num of indx_num before and after flipping
            vector<float> indx_num_af(indx_num);
            for (u32 i = 0; i < indx_be.size(); i++) {
                if ((u32)indx_be[i] != -1) {
                    indx_num_af[(u32)indx_be[i]] =
                        indx_num_af[(u32)indx_be[i]] - 1;
                }
            }
            for (u32 i = 0; i < indx_af.size(); i++) {
                if ((u32)indx_af[i] != -1) {
                    indx_num_af[(u32)indx_af[i]] =
                        indx_num_af[(u32)indx_af[i]] + 1;
                }
            }

            vector<float> change_indx;
            get_change_indx(indx_be, indx_af, change_indx);

            // Normalization between features and samples (single)
            PCA_mapper(indx_num_query, indx_num, indx_num_af, change_indx,
                       coeff, max_v, row_sum);

            // get the distance after flipping (single)
            vector<float> distance_af =
                get_distance(indx_num_af, w, norm_vector);

            // calculate the difference between distance_be and distance_af
            vector<float> difference(distance_be.size());
            for (u32 i = 0; i < distance_be.size(); i++) {
                difference[i] = abs(distance_be[i] - distance_af[i]);
            }

            // get the max difference
            float max_distance = (float)-9999999999;
            for (u32 i = 0; i < distance_be.size(); i++) {
                if (difference[i] > max_distance) {
                    max_distance = difference[i];
                }
            }
            dmap[i - edge_height][j - edge_width] = max_distance;
        }
    }
}

// int main(int argc, char const* argv[]) {
//     ifstream infile;
//     infile.open("img.txt");
//     vector<vector<float>> img(256, vector<float>(256, 0));
//     int k = 0;
//     for (u32 i = 0; i < img.size(); i++) {
//         for (u32 j = 0; j < img[0].size(); j++) {
//             infile >> img[i][j];
//         }
//     }
//     infile.close();

//     vector<u32> tmp_kern = {1,   2,   4,    8,    16,   32,   64,    128,
//                             256, 512, 1024, 2048, 4096, 8192, 16384, 32768};
//     vector<vector<float>> kern(4, vector<float>(4, 0));
//     k = 0;
//     for (u32 i = 0; i < kern.size(); i++) {
//         for (u32 j = 0; j < kern[i].size(); j++) {
//             kern[i][j] = tmp_kern[k++];
//         }
//     }

//     u32 dim_all = 5053;
//     u32 dim = 600;

//     infile.open("mu.txt");
//     vector<float> mu(dim_all);
//     k = 0;
//     for (u32 i = 0; i < mu.size(); i++) {
//         infile >> mu[i];
//     }
//     infile.close();

//     infile.open("coeff.txt");
//     vector<vector<float>> coeff(dim_all, vector<float>(dim, 0));
//     k = 0;
//     for (u32 i = 0; i < coeff.size(); i++) {
//         for (u32 j = 0; j < coeff[0].size(); j++) {
//             infile >> coeff[i][j];
//         }
//     }
//     infile.close();

//     infile.open("w.txt");
//     vector<vector<float>> w(1, vector<float>(dim, 0));
//     k = 0;
//     for (u32 i = 0; i < w.size(); i++) {
//         for (u32 j = 0; j < w[0].size(); j++) {
//             infile >> w[i][j];
//         }
//     }
//     infile.close();

//     infile.open("max_v.txt");
//     vector<float> max_v(dim_all);
//     k = 0;
//     for (u32 i = 0; i < max_v.size(); i++) {
//         infile >> max_v[i];
//     }
//     infile.close();

//     infile.open("pattern_select.txt");
//     vector<float> select(dim_all);
//     k = 0;
//     for (u32 i = 0; i < select.size(); i++) {
//         infile >> select[i];
//         select[i] -= 1;
//     }
//     infile.close();

//     infile.open("pattern_map.txt");
//     vector<float> pmap(65536);
//     k = 0;
//     for (u32 i = 0; i < pmap.size(); i++) {
//         infile >> pmap[i];
//         pmap[i] -= 1;
//     }
//     infile.close();

//     time_t start = clock();
//     vector<vector<float>> dmap(img.size(), vector<float>(img[0].size(), 0));
//     gen_distortion_map(dmap, img, kern, mu, coeff, w, max_v, select, pmap);

//     time_t stop = clock();

//     printf("%f seconds\n", (double)(stop - start) / CLOCKS_PER_SEC);

//     return 0;
// }

/** Matlat interface
GENDISTMAP calculation of distortion score map
        dmap = gendistmap(img, kern, mu, coeff, w, medv, disv)
---------------------------------------
IN
im - image vector<vector<float>>, float array
kern - kern filter
w - parameters of classification hyperplane
medv and disv - normalization between samples
OUT
dmap - distortion score map, float array
dmap_pair - distortion score map, float array
---------------------------------------
*/

void mexFunction(u32 nlhs, mxArray* plhs[], u32 nrhs, const mxArray* prhs[]) {

    if (nrhs != 8) {
        mexErrMsgTxt("The number of input parameters must 8.");
        return;
    }

    // IN img
    if (!mxIsSingle(prhs[MATLAB_IN_IMG])) {
        mexErrMsgTxt("The input image must be of type Single (float).");
        return;
    }
    vector<vector<float>> img(
        (u32)mxGetM(prhs[MATLAB_IN_IMG]),
        vector<float>((u32)mxGetN(prhs[MATLAB_IN_IMG]), 0));
    float* imptr = (float*)mxGetPr(prhs[MATLAB_IN_IMG]);
    for (u32 n = 0; n < img[0].size(); n++) {
        for (u32 m = 0; m < img.size(); m++) {
            img[m][n] = imptr[m + img.size() * n];
        }
    }

    // IN kern
    if (!mxIsSingle(prhs[MATLAB_IN_KERN])) {
        mexErrMsgTxt("The input kern must be of type Single (float).");
        return;
    }
    vector<vector<float>> kern(
        (u32)mxGetM(prhs[MATLAB_IN_KERN]),
        vector<float>((u32)mxGetN(prhs[MATLAB_IN_KERN]), 0));
    imptr = (float*)mxGetPr(prhs[MATLAB_IN_KERN]);
    for (u32 n = 0; n < kern[0].size(); n++) {
        for (u32 m = 0; m < kern.size(); m++) {
            kern[m][n] = imptr[m + kern.size() * n];
        }
    }

    // IN mu
    if (!mxIsSingle(prhs[MATLAB_IN_MU])) {
        mexErrMsgTxt("The input mu must be of type Single (float).");
        return;
    }
    if ((u32)mxGetM(prhs[MATLAB_IN_MU]) != 1) {
        mexErrMsgTxt("The input mu must be of type 1 x n.");
        return;
    }
    vector<float> mu((u32)mxGetN(prhs[MATLAB_IN_MU]), 0);
    imptr = (float*)mxGetPr(prhs[MATLAB_IN_MU]);
    for (u32 m = 0; m < mu.size(); m++) {
        mu[m] = imptr[m];
    }

    // IN coeff
    if (!mxIsSingle(prhs[MATLAB_IN_COEFF])) {
        mexErrMsgTxt("The input coeff must be of type Single (float).");
        return;
    }

    u32 select_dimension = (u32)mxGetM(prhs[MATLAB_IN_COEFF]);
    u32 last_dimension = (u32)mxGetN(prhs[MATLAB_IN_COEFF]);

    vector<vector<float>> coeff(select_dimension,
                                vector<float>(last_dimension, 0));
    imptr = (float*)mxGetPr(prhs[MATLAB_IN_COEFF]);
    for (u32 n = 0; n < coeff[0].size(); n++) {
        for (u32 m = 0; m < coeff.size(); m++) {
            coeff[m][n] = imptr[m + coeff.size() * n];
        }
    }

    // IN w
    if (!mxIsSingle(prhs[MATLAB_IN_W])) {
        mexErrMsgTxt("The input w must be of type Single (float).");
        return;
    }
    if ((u32)mxGetN(prhs[MATLAB_IN_W]) != last_dimension) {
        mexErrMsgTxt(
            "The input second dimension of w should be the same with "
            "pattern_select.");
        return;
    }

    u32 num_feature_space = (u32)mxGetM(prhs[MATLAB_IN_W]);

    vector<vector<float>> w(num_feature_space,
                            vector<float>(last_dimension, 0));
    imptr = (float*)mxGetPr(prhs[MATLAB_IN_W]);
    for (u32 n = 0; n < w[0].size(); n++) {
        for (u32 m = 0; m < w.size(); m++) {
            w[m][n] = imptr[m + w.size() * n];
        }
    }

    // IN max_v
    if (!mxIsSingle(prhs[MATLAB_IN_MAX_V])) {
        mexErrMsgTxt("The input max_v must be of type Single (float).");
        return;
    }
    if ((u32)mxGetN(prhs[MATLAB_IN_MAX_V]) != select_dimension) {
        mexErrMsgTxt(
            "The input second dimension of max_v should be the same "
            "with coeff.");
        return;
    }
    vector<float> max_v(select_dimension);
    imptr = (float*)mxGetPr(prhs[MATLAB_IN_MAX_V]);
    for (u32 m = 0; m < max_v.size(); m++) {
        max_v[m] = imptr[m];
    }

    // IN SELECT
    if (!mxIsSingle(prhs[MATLAB_IN_SELECT])) {
        mexErrMsgTxt(
            "The input pattern_select must be of type Single (float).");
        return;
    }
    if ((u32)mxGetM(prhs[MATLAB_IN_SELECT]) != 1) {
        mexErrMsgTxt("The input pattern_select must be of type 1 x n.");
        return;
    }
    if ((u32)mxGetN(prhs[MATLAB_IN_SELECT]) != select_dimension) {
        mexErrMsgTxt(
            "The input second dimension of pattern_select should be the same "
            "with coeff.");
        return;
    }
    vector<float> select(select_dimension);
    imptr = (float*)mxGetPr(prhs[MATLAB_IN_SELECT]);
    for (u32 m = 0; m < select.size(); m++) {
        select[m] = imptr[m] - 1;
    }

    // IN PMAP
    if (!mxIsSingle(prhs[MATLAB_IN_MAPPER])) {
        mexErrMsgTxt("The input pattern_map must be of type Single (float).");
        return;
    }
    if ((u32)mxGetN(prhs[MATLAB_IN_MAPPER]) != 65536) {
        mexErrMsgTxt("The input second dimension of pattern_map should 65536.");
        return;
    }
    if ((u32)mxGetM(prhs[MATLAB_IN_MAPPER]) != 1) {
        mexErrMsgTxt("The input pattern_map must be of type 1 x n.");
        return;
    }
    vector<float> pmap(65536);
    imptr = (float*)mxGetPr(prhs[MATLAB_IN_MAPPER]);
    for (u32 m = 0; m < pmap.size(); m++) {
        pmap[m] = imptr[m] - 1;
    }
    // end of MATLAB interface

    vector<vector<float>> dmap(img.size(), vector<float>(img[0].size(), 0));

    gen_distortion_map(dmap, img, kern, mu, coeff, w, max_v, select, pmap);
    mxArray* mDmap = mxCreateNumericMatrix(dmap.size(), dmap[0].size(),
                                           mxSINGLE_CLASS, mxREAL);
    imptr = (float*)mxGetPr(mDmap);
    for (u32 n = 0; n < dmap[0].size(); n++)
        for (u32 m = 0; m < dmap.size(); m++)
            imptr[m + dmap.size() * n] = dmap[m][n];

    plhs[MATLAB_OUT_DMAP] = mDmap;
}

#include "SImage.h"
#include "SImageIO.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <math.h>
#include <algorithm>
#include <limits.h>
#include <string>

using namespace std;

int DLIMIT = 255;

typedef enum {
	LEFT, RIGHT, UP, DOWN, ALL
} direction;

struct disp_arg {
	int start;
	int end;
	vector<SDoublePlane> &D;
	SDoublePlane &V;
	vector<vector<SDoublePlane> > &m;
	disp_arg(int a, int b, vector<SDoublePlane> &p, SDoublePlane &q,
			vector<vector<SDoublePlane> > &r) :
			start(a), end(b), D(p), V(q), m(r) {
	}
	;
};

// The simple image class is called SDoublePlane, with each pixel represented as
// a double (floating point) type. This means that an SDoublePlane can represent
// values outside the range 0-255, and thus can represent gradient magnitudes,
// harris corner scores, etc. 
//
// The SImageIO class supports reading and writing PNG files. It will read in
// a color PNG file, convert it to grayscale, and then return it to you in 
// an SDoublePlane. The values in this SDoublePlane will be in the range [0,255].
//
// To write out an image, call write_png_file(). It takes three separate planes,
// one for each primary color (red, green, blue). To write a grayscale image,
// just pass the same SDoublePlane for all 3 planes. In order to get sensible
// results, the values in the SDoublePlane should be in the range [0,255].
//

double get_messages_from_neighbors(SDoublePlane &m, int row, int col,
		direction dir) {
	double sum = 0;

	if (dir != LEFT)
		sum += m[row][col - 1];

	if (dir != RIGHT)
		sum += m[row][col + 1];

	if (dir != UP)
		sum += m[row - 1][col];

	if (dir != DOWN)
		sum += m[row + 1][col];

	return sum;
}

/*void set_disparity(void *obj)
 {
 disp_arg *A = (disp_arg*)obj;

 double disp_score;
 double min_score = INT_MAX;
 int label, t=1;
 double neighbors_sum;

 for (int i = A->start; i < A->end; i++) {
 for (int j = 1; j < A->D[0].cols() - 1; j++) {
 min_score = INT_MAX;
 label = DLIMIT;
 for (int d = 0; d < DLIMIT; d++) {
 neighbors_sum = get_messages_from_neighbors(A->m[d][t], i, j, ALL);
 disp_score = A->D[d][i][j] + neighbors_sum;

 if (disp_score < min_score) {
 min_score = disp_score;
 label = d;
 }
 }

 A->V[i][j] = label;
 }
 }
 }*/

void set_disparity(int start, int end, vector<SDoublePlane> &D, SDoublePlane &V,
		vector<vector<SDoublePlane> > &m) {
	//disp_arg *A = (disp_arg*)obj;

	double disp_score;
	double min_score = INT_MAX;
	int label, t = 1;
	double neighbors_sum;
	start = (start < 1) ? 1 : start;
	end = (end < (V.rows() - 1)) ? end : (V.rows() - 1);
	for (int i = start; i < end - 1; i++) {
		for (int j = 1; j < D[0].cols() - 1; j++) {
			min_score = INT_MAX;
			label = DLIMIT;
			for (int d = 0; d < DLIMIT; d++) {
				//cout << i<<","<<j<<endl;
				neighbors_sum = get_messages_from_neighbors(m[d][t], i, j, ALL);
				disp_score = D[d][i][j] + neighbors_sum;

				if (disp_score < min_score) {
					min_score = disp_score;
					label = d;
				}
			}

			V[i][j] = label;
		}
	}
}

void compute_unary_cost(vector<SDoublePlane> &D, const SDoublePlane &left_image,
		const SDoublePlane &right_image, SDoublePlane &disp, int w) {
	double sum = 0;
	for (int i = w; i < left_image.rows() - w; i++) {
		for (int j = w; j < left_image.cols() - w; j++) {
			for (int d = 0; d < DLIMIT; d++) {
				sum = 0;
				for (int u = -w; u < w; u++) {
					for (int v = -w; v < w; v++) {
						int idx =
								((j + v + d) > left_image.cols()) ?
										left_image.cols() : j + v + d;
						sum += pow(
								left_image[i + u][j + v]
										- right_image[i + u][idx], 2);
					}
				}
				// find label based on minimum difference
//			disp[i][j] = D[disp[i][j]][i][j] <= sum ? disp[i][j] : d;
//			cout<<disp[i][j]<<endl;
				D[d][i][j] = sum / w * w;
			}

		}
	}

}

double compute_pairwise_cost(int i, int j, SDoublePlane &disp) {
	int sum = 0;
	int min_diff = INT_MAX;
	int min_label = DLIMIT;

//	SDoublePlane result(disp.rows(), disp.cols());

//	for (int i = 1; i < disp.rows() - 1; i++) {
//		for (int j = 1; j < disp.cols() - 1; j++) {
	min_diff = INT_MAX;

	for (int d = 0; d < DLIMIT; d++) {
		sum = 0;
		sum += pow(d - disp[i][j - 1], 2);
		sum += pow(d - disp[i][j + 1], 2);
		sum += pow(d - disp[i + 1][j], 2);
		sum += pow(d - disp[i - 1][j], 2);
		if (sum < min_diff) {
			min_diff = sum;
			min_label = d;
		}
	}
	return min_diff;
//			result[i][j] = min_diff;
//			disp[i][j] = min_label;
	//cout<<min_diff<<endl;

}

void propogate_belief(vector<SDoublePlane> &D, SDoublePlane &V,
		const SDoublePlane &left_image, const SDoublePlane &right_image) {

	vector<vector<SDoublePlane> > m;
	vector<SDoublePlane> tmp;

	SDoublePlane probs(D[0].rows(), D[0].cols());

	tmp.push_back(probs);
	tmp.push_back(probs);

	for (int i = 0; i < DLIMIT; i++)
		m.push_back(tmp);

	int t = 1;
	int t_minus_one = 0;
	int target_label = 0;
	int source_label = DLIMIT;
	double no_match_score = INT_MAX;
	double neighbors_sum = 0;
	int C = 1;
	//double src_nbrs_sum = 0;

	cout << "Starting loopy bp \n";
	int it = 0;
	compute_unary_cost(D, left_image, right_image, V, 3);
	while (it++ < 1) {

		cout << "iterations " << it << endl;
		for (int i = 1; i < probs.rows() - 1; i++) {
			for (int j = probs.cols() - 1; j > 0; j--) {
				double min_score = INT_MAX;


				for (source_label = 0; source_label < DLIMIT - 1; source_label++) {
					neighbors_sum = get_messages_from_neighbors(
							m[source_label][t_minus_one], i, j, LEFT);
					min_score = min(min_score,
							neighbors_sum + D[source_label][i][j]);
				}

				for (target_label = 0; target_label < DLIMIT - 1;
						target_label++) {
					if (target_label == source_label) {
						neighbors_sum = get_messages_from_neighbors(
								m[source_label][t_minus_one], i, j, LEFT);
						m[target_label][t][i][j - 1] = neighbors_sum
								+ D[source_label][i][j];
					} else {
						m[target_label][t][i][j - 1] = min_score + C;
					}
				}
			}
		}
		//		cout<<m[target_label].size()<<endl;
		cout << "Messages sent left\n";

		for (int i = 1; i < probs.rows() - 1; i++) {
			for (int j = 1; j < probs.cols() - 1; j++) {
				double min_score = INT_MAX;

				for (source_label = 0; source_label < DLIMIT - 1; source_label++) {
					neighbors_sum = get_messages_from_neighbors(
							m[source_label][t_minus_one], i, j, RIGHT);
					min_score = min(min_score,
							neighbors_sum + D[source_label][i][j]);
				}

				for (target_label = 0; target_label < DLIMIT - 1;
						target_label++) {
					if (target_label == source_label) {
						neighbors_sum = get_messages_from_neighbors(
								m[source_label][t_minus_one], i, j, RIGHT);
						m[target_label][t][i][j + 1] = neighbors_sum
								+ D[source_label][i][j];
					} else {
						m[target_label][t][i][j + 1] = min_score + C;
					}
				}
			}
		}

		cout << "Messages sent right\n";

		int tmp = t;
		t = t_minus_one;
		t_minus_one = t;

	}

	double disp_score;
	double min_score = INT_MAX;
	int label;

	for (int i = 1; i < probs.rows() - 1; i++) {
		for (int j = 1; j < probs.cols() - 1; j++) {
			min_score = INT_MAX;
			label = DLIMIT;
			for (int d = 0; d < DLIMIT; d++) {
				neighbors_sum = get_messages_from_neighbors(m[d][t], i, j, ALL);
				disp_score = D[d][i][j] + neighbors_sum;
				if (disp_score < min_score) {
					min_score = disp_score;
					label = d;
				}
			}
			V[i][j] = 5* label;
		}
	}
}

SDoublePlane mrf_stereo(const SDoublePlane &left_image,
		const SDoublePlane &right_image) {
	// implement this in step 4...
	//  this placeholder just returns a random disparity map
	SDoublePlane result(left_image.rows(), left_image.cols());

	// disparity level for each pixel
	SDoublePlane disp(left_image.rows(), left_image.cols());
	SDoublePlane labels(left_image.rows(), left_image.cols());
	vector<SDoublePlane> D;

	DLIMIT = 50;// left_image.cols() - 1;

	for (int i = 0; i < left_image.rows(); i++) {
		for (int j = 0; j < left_image.cols(); j++) {
			disp[i][j] = DLIMIT - 1;
		}
	}

	for (int i = 0; i < left_image.rows(); i++) {
		for (int j = 0; j < left_image.cols(); j++) {
			labels[i][j] = INT_MAX;
		}
	}

	for (int i = 0; i < DLIMIT; i++)
		D.push_back(labels);

	for (int iter = 0; iter < 1; iter++) {
//		for (int i = 0; i < DLIMIT; i++) {
//			compute_unary_cost(D, left_image, right_image, disp, 5);
//		}
//		result = compute_pairwise_cost(disp);
		propogate_belief(D, disp, left_image, right_image);
	}

	/*for(int i=0; i<left_image.rows(); i++)
	 for(int j=0; j<left_image.cols(); j++)
	 result[i][j] = rand() % 256;*/

	return disp;
}

int main(int argc, char *argv[]) {
	if (argc != 4 && argc != 3) {
		cerr << "usage: " << argv[0] << " image_file1 image_file2 [gt_file]"
				<< endl;
		return 1;
	}

	string input_filename1 = argv[1], input_filename2 = argv[2];
	string gt_filename;
	if (argc == 4)
		gt_filename = argv[3];

	// read in images and gt
	SDoublePlane image1 = SImageIO::read_png_file(input_filename1.c_str());
	SDoublePlane image2 = SImageIO::read_png_file(input_filename2.c_str());
	SDoublePlane gt;
	if (gt_filename != "") {
		gt = SImageIO::read_png_file(gt_filename.c_str());
		// gt maps are scaled by a factor of 3, undo this...
		for (int i = 0; i < gt.rows(); i++)
			for (int j = 0; j < gt.cols(); j++)
				gt[i][j] = gt[i][j] / 3.0;
	}

	// do stereo using mrf
	SDoublePlane disp3 = mrf_stereo(image1, image2);
	SImageIO::write_png_file("disp_mrf.png", disp3, disp3, disp3);

	// Measure error with respect to ground truth, if we have it...
	if (gt_filename != "") {
		double err = 0;
		for (int i = 0; i < gt.rows(); i++)
			for (int j = 0; j < gt.cols(); j++)
				err += sqrt(
						(disp3[i][j] - gt[i][j]) * (disp3[i][j] - gt[i][j]));

		cout << "MRF stereo technique mean error = "
				<< err / gt.rows() / gt.cols() << endl;

	}

	return 0;
}

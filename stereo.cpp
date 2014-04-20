#include "SImage.h"
#include "SImageIO.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <math.h>
#include <algorithm>
#include <limits.h>
#include<string>

#define WINDOW_SIZE 2


using namespace std;

int DLIMIT = 100;

typedef enum {
	LEFT, RIGHT, UP, DOWN, ALL
} direction;



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


void display(const SDoublePlane &d)
{
	for(int i=0; i<d.rows(); i++)
	{
		for(int j=0; j<d.cols(); j++)
		{
			cout<<d[i][j]<<",";
		}
		cout<<endl;
	}
}



double get_messages_from_neighbors(SDoublePlane &m, int row, int col,
		direction dir) {
	double sum = 0;
	int n = 0;

	if (dir != LEFT) {
		sum += m[row][col - 1];
		n++;
	}

	if (dir != RIGHT) {
		sum += m[row][col + 1];
		n++;
	}

	if (dir != UP) {
		sum += m[row - 1][col];
		n++;
	}

	if (dir != DOWN) {
		sum += m[row + 1][col];
		n++;
	}

	return (sum/n);
}

/*void compute_unary_cost(const SDoublePlane &left_image,
		const SDoublePlane &right_image, vector<SDoublePlane> &D,SDoublePlane &disp) {
	double sum = 0;
	int w = WINDOW_SIZE;

	for (int i = w; i < left_image.rows() - w; i++) {
		for (int j = w; j < left_image.cols() - w; j++) {
			for (int d = 0; d < DLIMIT; d++) {
			sum = 0;
			for (int u = -w; u < w; u++) {
				for (int v = -w; v < w; v++) {
					int idx = ((j+v+d) > left_image.cols())? left_image.cols():j+v+d;
					sum += pow(	left_image[i + u][j + v] - right_image[i + u][idx], 2);
				}
			}
			// find label based on minimum difference
			disp[i][j] = D[disp[i][j]][i][j] <= sum ? disp[i][j] : d;

			D[d][i][j] = sum/w*w;
			}


		}
	}

}*/


void compute_unary_cost(const SDoublePlane &left_image, const SDoublePlane &right_image,
			vector<SDoublePlane> &D, SDoublePlane &V) {
	double sum = 0;
	int w = WINDOW_SIZE;
	int val = 0;
	int row = 0;
	int col = 0;

	cout<<"Computing Unary Potential..."<<endl;

	for (int d = 0; d < DLIMIT; d++) {
		for (int i = w; i <= left_image.rows() - w - 1; i++) {
			for (int j = w; j <= left_image.cols() - w - 1; j++) {

				sum = 0;
				val = 0;

				for (int u = -w; u <= w; u++) {
					for (int v = -w; v <= w; v++) {

						row = i + u;
						col = j + v;

						if (col + d < right_image.cols() - 1)
							val = right_image[row][col + d];

						sum += pow(left_image[row][col] - val, 2);
					}
				}

				sum = sum /(4*w*w);

				if (sum <= D[V[i][j]][i][j]) {
					V[i][j] = d;
				}
				D[d][i][j] = sum;
			}
		}
	}
}





void propogate_belief(vector<SDoublePlane> &D, SDoublePlane &V,
		const SDoublePlane &left_image,
		const SDoublePlane &right_image) {

	vector<vector<SDoublePlane> >m;
	vector<SDoublePlane> tmp;
	//vector<vector<SDoublePlane> > dtmp;


	SDoublePlane probs(D[0].rows(), D[0].cols());

	for(int i=0; i<probs.rows(); i++)
		for(int j=0; j<probs.cols(); j++)
			probs[i][j] = 0;

	tmp.push_back(probs);
	tmp.push_back(probs);

	for(int i=0; i<DLIMIT; i++)
		m.push_back(tmp);

	int t = 1;
	int t_minus_one = 0;
	int target_label = 0;
	int source_label = DLIMIT;
	int potts_cost = DLIMIT/2;

	double tmp_score = 0;
	double match_score = 0;
	double no_match_score = 0;
	double neighbors_sum = 0;
	double disp_score;
	double min_score = INT_MAX;


	int it = 0;

	while (it++ < 1) {

		cout<<"Iteration : "<<it<<"\n------------------\n";

		for (int i = 1; i < probs.rows() - 1; i++) {
			for (int j = probs.cols() - 1; j > 0; j--) {

				min_score = INT_MAX;
				source_label = 0;

				for(int qi=0; qi<DLIMIT; qi++)
				{
					neighbors_sum = get_messages_from_neighbors(m[qi][t_minus_one], i, j, LEFT);
					tmp_score = D[qi][i][j] + neighbors_sum;

					if(tmp_score < min_score)
					{
						min_score = tmp_score;
						source_label = qi;
					}
				}

				///no_match_score = INT_MAX;
				//match_score = INT_MAX;
				//target_label = 0;

				for(int qj=0; qj<DLIMIT; qj++)
				{
					if(qj == source_label) continue;

					neighbors_sum = get_messages_from_neighbors(m[qj][t_minus_one], i, j, LEFT);
					match_score = D[qj][i][j] + neighbors_sum;

					m[qj][t][i][j-1] = min(match_score, min_score + pow(qj-source_label,2));

				}


				//m[target_label][t][i][j-1] = min(match_score, no_match_score);
			}
		}

		cout << "Messages sent left...\n";

		for (int i = 1; i < probs.rows() - 1; i++) {
			for (int j = 1; j < probs.cols() - 1; j++) {

				min_score = INT_MAX;
				source_label = 0;

				for (int qi = 0; qi < DLIMIT; qi++) {
					neighbors_sum = get_messages_from_neighbors(
							m[qi][t_minus_one], i, j, RIGHT);
					tmp_score = D[qi][i][j] + neighbors_sum;

					if (tmp_score < min_score) {
						min_score = tmp_score;
						source_label = qi;
					}
				}

				for (int qj = 0; qj < DLIMIT; qj++) {

					if(qj == source_label) continue;

					neighbors_sum = get_messages_from_neighbors(
							m[qj][t_minus_one], i, j, RIGHT);
					match_score = D[qj][i][j] + neighbors_sum;

					m[qj][t][i][j - 1] = min(match_score,
							min_score + pow(qj - source_label, 2));

				}

				//m[target_label][t][i][j+1] = min(match_score, no_match_score);
			}
		}

		cout << "Messages sent right...\n";


		for (int j = 1; j < probs.cols() - 1; j++) {
			for (int i = probs.rows() - 2; i > 0; i--) {

				min_score = INT_MAX;
				source_label = 0;

				for (int qi = 0; qi < DLIMIT; qi++) {
					neighbors_sum = get_messages_from_neighbors(
							m[qi][t_minus_one], i, j, UP);
					tmp_score = D[qi][i][j] + neighbors_sum;

					if (tmp_score < min_score) {
						min_score = tmp_score;
						source_label = qi;
					}
				}

				for (int qj = 0; qj < DLIMIT; qj++) {

					if (qj == source_label)
						continue;

					neighbors_sum = get_messages_from_neighbors(
							m[qj][t_minus_one], i, j, UP);
					match_score = D[qj][i][j] + neighbors_sum;

					m[qj][t][i-1][j] = min(match_score,
							min_score + pow(qj - source_label, 2));

				}
			}
		}

		cout << "Messages sent up...\n";

		for (int j = 1; j < probs.cols() - 1; j++) {
			for (int i = 1; i < probs.rows() - 1; i++) {

				min_score = INT_MAX;
				source_label = 0;

				for (int qi = 0; qi < DLIMIT; qi++) {
					neighbors_sum = get_messages_from_neighbors(
							m[qi][t_minus_one], i, j, DOWN);
					tmp_score = D[qi][i][j] + neighbors_sum;

					if (tmp_score < min_score) {
						min_score = tmp_score;
						source_label = qi;
					}
				}

				for (int qj = 0; qj < DLIMIT; qj++) {

					if (qj == source_label)
						continue;

					neighbors_sum = get_messages_from_neighbors(
							m[qj][t_minus_one], i, j, DOWN);
					match_score = D[qj][i][j] + neighbors_sum;

					m[qj][t][i+1][j] = min(match_score,
							min_score + pow(qj - source_label, 2));

				}
			}
		}
		cout << "Messages sent down...\n";

		cout<<"------------------\n";


		int tmp = t;
		t = t_minus_one;
		t_minus_one = t;
	}


	cout << "Updating Labels...\n";

	for (int i = 1; i < probs.rows() - 1; i++) {
		for (int j = 1; j < probs.cols() - 1; j++) {
			min_score = INT_MAX;

			for (int d = 0; d < DLIMIT; d++) {
				neighbors_sum = get_messages_from_neighbors(m[d][t], i, j, ALL);

				disp_score = D[d][i][j] + neighbors_sum;

				if (disp_score < min_score) {
					min_score = disp_score;
					V[i][j] = d;
				}
			}
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



	for (int i = 0; i < left_image.rows(); i++) {
		for (int j = 0; j < left_image.cols(); j++) {
			disp[i][j] = INT_MAX;
			labels[i][j] = DLIMIT - 1;
		}
	}

	for (int i = 0; i < DLIMIT; i++)
		D.push_back(disp);

	cout << "MAX DISPARITY : " << DLIMIT << endl;
	cout << "WINDOW SIZE : " << WINDOW_SIZE << endl;

	compute_unary_cost(left_image, right_image, D, labels);

	//display(labels);

	propogate_belief(D, labels, left_image, right_image);

	for(int i=0; i<labels.rows(); i++)
		for(int j=0; j<labels.cols(); j++)
			labels[i][j] = labels[i][j] * 256/DLIMIT;

	return labels;
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

	cout<<"--Done--\n";

	return 0;
}

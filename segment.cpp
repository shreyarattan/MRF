#include "SImage.h"
#include "SImageIO.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <math.h>
#include <algorithm>
#include <limits.h>

#define FG 1
#define BG 0

typedef enum {
	LEFT, RIGHT, UP, DOWN, ALL
} direction;

using namespace std;

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

class Point {
public:
	Point() {
	}
	Point(int _r, int _c) :
			row(_r), col(_c) {
	}
	int row, col;
};

bool find_point(const vector<Point> fg, Point p) {
	for (vector<Point>::const_iterator iter = fg.begin(); iter != fg.end();
			++iter)
		if (iter->col == p.col && iter->row == p.row) {
			return true;
		}
	return false;
}

void calc_mean(SDoublePlane red_plane, SDoublePlane green_plane,
		SDoublePlane blue_plane, vector<double>&m, vector<double>&sig,
		const vector<Point> &fg) {
	double sum_red = 0, sum_green = 0, sum_blue = 0;
	for (int i = 0; i < red_plane.rows(); i++)
		for (int j = 0; j < red_plane.cols(); j++) {
			Point p(i, j);
			if (find_point(fg, p)) {
				sum_red += red_plane[i][j];
				sum_green += green_plane[i][j];
				sum_blue += blue_plane[i][j];
			}
		}

	double mean_red = sum_red / (fg.size());
	double mean_green = sum_green / (fg.size());
	double mean_blue = sum_blue / (fg.size());
	m.push_back(mean_red);
	m.push_back(mean_green);
	m.push_back(mean_blue);

	sum_red = 0, sum_green = 0, sum_blue = 0;

	for (int i = 0; i < red_plane.rows(); i++)
		for (int j = 0; j < red_plane.cols(); j++) {
			Point p(i, j);
			if (find_point(fg, p)) {
				sum_red += pow((red_plane[i][j] - mean_red), 2);
				sum_green += pow((green_plane[i][j] - mean_green), 2);
				sum_blue += pow((blue_plane[i][j] - mean_blue), 2);
			}
		}

	double variance = sum_red / (fg.size() - 1);
	sig.push_back(variance);
	variance = sum_green / (fg.size() - 1);
	sig.push_back(variance);
	variance = sum_blue / (fg.size() - 1);
	sig.push_back(variance);

}

double get_gaussian_probs(double red_val, double green_val, double blue_val,
		vector<double> m, vector<double> sig) {
	return (exp(-(pow(red_val - m[0], 2) / (pow(sig[0], 2))))
			* exp(-(pow(green_val - m[1], 2) / (pow(sig[1], 2))))
			* exp(-(pow(blue_val - m[2], 2) / (pow(sig[2], 2)))));
}

double get_gaussian(double red_val, double green_val, double blue_val,
		vector<double> m, vector<double> sig) {
	return (-(pow(red_val - m[0], 2) / (pow(sig[0], 2)))
			* -(pow(green_val - m[1], 2) / (pow(sig[1], 2)))
			* -(pow(blue_val - m[2], 2) / (pow(sig[2], 2))));
}

SDoublePlane probs_table;

SDoublePlane naive_segment(const SDoublePlane *img, const vector<Point> &fg,
		const vector<Point> &bg) {
	// implement this in step 2...
	//  this placeholder just returns a silly segmentation (splitting the image in 2 arbitrarily)

	// To make things a little easier for you, we'll create separate red, green, and blue
	//  planes from the input image.

	const SDoublePlane &red_plane = img[0], &green_plane = img[1], blue_plane =
			img[2];
	vector<double> m, sig;

	calc_mean(red_plane, green_plane, blue_plane, m, sig, fg);

	SDoublePlane result(red_plane.rows(), red_plane.cols());
	double fg_score = INT_MAX;
	double min_probability = INT_MAX;
	double max_probability = INT_MIN;

	probs_table = SDoublePlane(red_plane.rows(), red_plane.cols());

	double beta = 0.97;

	for (int i = 0; i < result.rows(); i++) {
		for (int j = 0; j < result.cols(); j++) {
			Point p(i, j);
			if (find_point(fg, p)) {
				result[i][j] = 1;
			} else if (find_point(bg, p)) {
				result[i][j] = 0;
			} else if (!find_point(fg, p) && !find_point(bg, p)) {

				fg_score = get_gaussian_probs(red_plane[i][j],
						green_plane[i][j], blue_plane[i][j], m, sig);

				if (fg_score > beta)
					result[i][j] = 1;
				else
					result[i][j] = 0;

				probs_table[i][j] = (fg_score > beta) ? fg_score : beta;

				if (fg_score > max_probability)
					max_probability = fg_score;
				if (fg_score < min_probability)
					min_probability = fg_score;

			}
		}
	}

	cout << "MAX :" << max_probability << ", MIN :" << min_probability << endl;

	return result;
}

double get_messages_from_neighbors(SDoublePlane m, int row, int col,
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

void propogate_belief(vector<SDoublePlane> D, SDoublePlane &V) {
	vector<vector<SDoublePlane> > m;
	vector<SDoublePlane> tmp;

	SDoublePlane probs(D[FG].rows(), D[FG].cols());

	tmp.push_back(probs);
	tmp.push_back(probs);

	m.push_back(tmp);
	m.push_back(tmp);

	int t = 1;
	int t_minus_one = 0;
	int target_label = BG;
	int source_label = FG;
	int potts_cost = 1;

	double FGScore = 0;
	double BGScore = 0;
	double match_score = 0;
	double no_match_score = 0;
	double neighbors_sum = 0;
	//double src_nbrs_sum = 0;
	double tmp_diff = 0;
	double diff = 100;

	cout << "Starting loopee bee pee \n";

	int it = 0;

	while (it++ < 3) {
		for (int i = 1; i < probs.rows() - 1; i++) {
			for (int j = 1; j < probs.cols() - 1; j++) {

				target_label = V[i][j + 1];
				source_label = V[i][j];

				neighbors_sum = get_messages_from_neighbors(
						m[target_label][t_minus_one], i, j, RIGHT);
				match_score = D[target_label][i][j] + neighbors_sum;

				neighbors_sum = get_messages_from_neighbors(
						m[source_label][t_minus_one], i, j, RIGHT);
				no_match_score = potts_cost + D[source_label][i][j]
						+ neighbors_sum;

				m[target_label][t][i][j + 1] = min(match_score, no_match_score);

				//if pixel was FG
				/*neighbors_sum = get_messages_from_neighbors(m[FG][t_minus_one], i, j, RIGHT);

				 FGScore = D[FG][i][j] +
				 pow(FG - V[i][j+1], 2) +
				 neighbors_sum;

				 neighbors_sum = get_messages_from_neighbors(m[BG][t_minus_one], i, j, RIGHT);

				 BGScore = D[BG][i][j] +
				 pow(BG - V[i][j+1], 2) +
				 neighbors_sum;

				 m[V[i][j+1]][t][i][j+1] = min(FGScore, BGScore);*/

				//tmp_diff += m[V[i][j+1]][t][i][j+1] - m[V[i][j+1]][t_minus_one][i][j+1];
			}
		}

		cout << "Messages sent right\n";

		//diff = tmp_diff/(probs.rows()*probs.cols());
		//tmp_diff = 0;

		for (int i = 1; i < probs.rows() - 1; i++) {
			for (int j = probs.cols() - 1; j > 0; j--) {

				/*target_label = V[i][j-1];
				 source_label = V[i][j];

				 neighbors_sum = get_messages_from_neighbors(m[target_label][t_minus_one], i, j, LEFT);
				 match_score = D[target_label][i][j] + neighbors_sum;

				 neighbors_sum = get_messages_from_neighbors(m[source_label][t_minus_one], i, j, LEFT);
				 no_match_score = potts_cost + D[source_label][i][j] + neighbors_sum;

				 m[target_label][t][i][j-1] = min(match_score, no_match_score);		*/

				//if pixel was FG
				neighbors_sum = get_messages_from_neighbors(m[FG][t_minus_one],
						i, j, LEFT);

				FGScore = D[FG][i][j] + pow(FG - V[i][j - 1], 2)
						+ neighbors_sum;

				neighbors_sum = get_messages_from_neighbors(m[BG][t_minus_one],
						i, j, LEFT);

				BGScore = D[BG][i][j] + pow(BG - V[i][j - 1], 2)
						+ neighbors_sum;

				m[V[i][j - 1]][t][i][j - 1] = min(FGScore, BGScore);

				//tmp_diff += m[V[i][j-1]][t][i][j+1] - m[V[i][j-1]][t_minus_one][i][j+1];
			}
		}

		cout << "Messages sent left\n";

		//tmp_diff = tmp_diff/(probs.rows()*probs.cols());
		//diff = (diff + tmp_diff)/2;

		cout << "Energy diff :" << diff << endl;

		int tmp = t;
		t = t_minus_one;
		t_minus_one = t;
	}

	for (int i = 1; i < probs.rows() - 1; i++) {
		for (int j = 1; j < probs.cols() - 1; j++) {

			neighbors_sum = get_messages_from_neighbors(m[FG][t], i, j, ALL);
			FGScore = D[FG][i][j] + neighbors_sum;

			neighbors_sum = get_messages_from_neighbors(m[BG][t], i, j, ALL);
			BGScore = D[BG][i][j] + neighbors_sum;

			V[i][j] = FGScore < BGScore ? FG : BG;
		}
	}

}

SDoublePlane mrf_segment(const SDoublePlane *img, const vector<Point> &fg,
		const vector<Point> &bg) {
	// implement this in step 3...
	//  this placeholder just returns a random disparity map by calling naive_segment

	SDoublePlane fg_energy(img->rows(), img->cols());
	SDoublePlane bg_energy(img->rows(), img->cols());
	SDoublePlane result(img->rows(), img->cols());
	vector<SDoublePlane> D;
	SDoublePlane V;

	const SDoublePlane &red_plane = img[0], &green_plane = img[1], blue_plane =
			img[2];
	vector<double> m, sig;
//
//    double beta = 2e-009;
//
	calc_mean(red_plane, green_plane, blue_plane, m, sig, fg);

	double fg_score = INT_MAX;
	double min_probability = INT_MAX;
	double max_probability = INT_MIN;

	probs_table = SDoublePlane(red_plane.rows(), red_plane.cols());

	double beta = 1e-008;

	for (int i = 0; i < result.rows(); i++) {
		for (int j = 0; j < result.cols(); j++) {
			Point p(i, j);
			if (find_point(fg, p)) {
				result[i][j] = 1;
				fg_energy[i][j] = 0;
				bg_energy[i][j] = INT_MAX;
			} else if (find_point(bg, p)) {
				result[i][j] = 0;
				bg_energy[i][j] = 0;
				fg_energy[i][j] = INT_MAX;
			} else if (!find_point(fg, p) && !find_point(bg, p)) {
				fg_energy[i][j] = -(get_gaussian(red_plane[i][j],
						green_plane[i][j], blue_plane[i][j], m, sig));
				bg_energy[i][j] = beta;

				result[i][j] = fg_energy[i][j] < bg_energy[i][j] ? FG : BG;
			}
		}
	}

	D.push_back(bg_energy);
	D.push_back(fg_energy);

	propogate_belief(D, result);

	//result = loopy_belief(v, probs_table);

	cout << "MRF Done" << endl;

	return result;
}

// Take in an input image and a binary segmentation map. Use the segmentation map to split the 
//  input image into foreground and background portions, and then save each one as a separate image.
//
// This code should work okay already, but feel free to modify if you want.
//
void output_segmentation(const SDoublePlane *img, const SDoublePlane &labels,
		const string &fname) {
	// sanity checks. If one of these asserts fails, you've given this function invalid arguments!
	assert(img[0].rows() == labels.rows());
	assert(img[0].cols() == labels.cols());

	SDoublePlane img_fg[3], img_bg[3];

	for (int i = 0; i < 3; i++) {
		img_fg[i] = img[i];
		img_bg[i] = img[i];
	}

	for (int i = 0; i < labels.rows(); i++)
		for (int j = 0; j < labels.cols(); j++) {
			if (labels[i][j] == 0)
				img_fg[0][i][j] = img_fg[1][i][j] = img_fg[2][i][j] = 0;
			else if (labels[i][j] == 1)
				img_bg[0][i][j] = img_bg[1][i][j] = img_bg[2][i][j] = 0;
			else
				assert(0);
		}

	SImageIO::write_png_file((fname + "_fg.png").c_str(), img_fg[0], img_fg[1],
			img_fg[2]);
	SImageIO::write_png_file((fname + "_bg.png").c_str(), img_bg[0], img_bg[1],
			img_bg[2]);
}

int main(int argc, char *argv[]) {
	if (argc != 3) {
		cerr << "usage: " << argv[0] << " image_file seeds_file" << endl;
		return 1;
	}

	string input_filename1 = argv[1], seeds_file = argv[2];

	// read in images and gt
	SDoublePlane image_rgb[3], seeds_rgb[3];
	SImageIO::read_png_file_rgb(input_filename1.c_str(), image_rgb);
	SImageIO::read_png_file_rgb(seeds_file.c_str(), seeds_rgb);

	// figure out seed points
	vector<Point> fg_pixels, bg_pixels;
	for (int i = 0; i < seeds_rgb[0].rows(); i++)
		for (int j = 0; j < seeds_rgb[0].cols(); j++) {
			// blue --> foreground
			if (seeds_rgb[0][i][j] < 100 && seeds_rgb[1][i][j] < 100
					&& seeds_rgb[2][i][j] > 100)
				fg_pixels.push_back(Point(i, j));

			// red --> background
			if (seeds_rgb[0][i][j] > 100 && seeds_rgb[1][i][j] < 100
					&& seeds_rgb[2][i][j] < 100)
				bg_pixels.push_back(Point(i, j));
		}

	// do naive segmentation

//	SDoublePlane labels = naive_segment(image_rgb, fg_pixels, bg_pixels);
//	output_segmentation(image_rgb, labels, "naive_segment_result");

	// do mrf segmentation
	SDoublePlane labels = mrf_segment(image_rgb, fg_pixels, bg_pixels);
	output_segmentation(image_rgb, labels, "mrf_segment_result");

	return 0;
}

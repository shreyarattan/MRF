#include <SImage.h>
#include <SImageIO.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <math.h>
#include <algorithm>

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
	Point() { }
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
			if (find_point(fg, p)){
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

double get_gaussian(double red_val, double green_val, double blue_val,
		vector<double> m, vector<double> sig) {
	return (exp(-(pow(red_val - m[0], 2) / (pow(sig[0], 2))))
			* exp(-(pow(green_val - m[1], 2) / (pow(sig[1], 2))))
			* exp(-(pow(blue_val - m[2], 2) / (pow(sig[2], 2)))));
}

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
//	cout << "red:" << m[0] << "," << sig[0] << endl;
//	cout << "green:" << m[1] << "," << sig[1] << endl;
//	cout << "blue:" << m[2] << "," << sig[2] << endl;

	SDoublePlane result_0(red_plane.rows(), red_plane.cols());
	SDoublePlane result_1(red_plane.rows(), red_plane.cols());

	for (int i = 0; i < red_plane.rows(); i++) {
		for (int j = 0; j < red_plane.cols(); j++) {
			result_0[i][j] = 0;
			result_1[i][j] = 1;
		}
	}

	//assign forground and background colors according to the colors strokes
	/*for(vector<Point>::const_iterator iter = fg.begin(); iter != fg.end(); ++iter)
	 result_0[iter->row][iter->col] = 1;

	 for(vector<Point>::const_iterator iter = bg.begin(); iter != bg.end(); ++iter)
	 result_1[iter->row][iter->col] = 0;*/
	double beta = 0.9999000;
	cout << beta<<endl;
	for (int i = 0; i < result_0.rows(); i++) {
		for (int j = 0; j < result_0.cols(); j++) {
			Point p(i, j);
			if (find_point(fg, p)) {
				result_0[i][j] = 0;
				result_1[i][j] = 9999999;
			} else if (find_point(bg, p)) {
				result_1[i][j] = 0;
				result_0[i][j] = 9999999;
			} else if (!find_point(fg, p) && !find_point(bg, p)) {
				if (result_0[i][j] == 0)
					result_0[i][j] = beta;
				if (result_1[i][j] == 1) {
					result_1[i][j] = get_gaussian(red_plane[i][j],
							green_plane[i][j], blue_plane[i][j], m, sig);
				}
			}
		}
	}
	SDoublePlane result(red_plane.rows(), red_plane.cols());
	for (int i = 0; i < result_0.rows(); i++) {
		for (int j = 0; j < result_0.cols(); j++) {
//		if (result_0[i][j] == 1) result[i][j] = 1;
//		else if (result_1[i][j] == 0) result[i][j] = 0;
//		else {
			result[i][j] = result_0[i][j] < result_1[i][j] ? 1 : 0;
			cout << result_0[i][j] << "," << result_1[i][j] << "," << result[i][j] << endl;
//		}
		}
	}

	// Now you can use red_plane, green_plane, and blue_plane just like any other planes.
	//SDoublePlane result(red_plane.rows(), red_plane.cols());

	/*for(vector<Point>::const_iterator iter = fg.begin(); iter != fg.end(); ++iter)
	 result[iter->row][iter->col] = 1;

	 for(vector<Point>::const_iterator iter = bg.begin(); iter != bg.end(); ++iter)
	 result[iter->row][iter->col] = 0;
	 */

	return result;
}

SDoublePlane mrf_segment(const SDoublePlane *img, const vector<Point> &fg,
		const vector<Point> &bg) {
	// implement this in step 3...
	//  this placeholder just returns a random disparity map by calling naive_segment
	return naive_segment(img, fg, bg);
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
	SDoublePlane labels = naive_segment(image_rgb, fg_pixels, bg_pixels);
	output_segmentation(image_rgb, labels, "naive_segment_result");

	// do mrf segmentation
	labels = naive_segment(image_rgb, fg_pixels, bg_pixels);
	output_segmentation(image_rgb, labels, "mrf_segment_result");

	return 0;
}

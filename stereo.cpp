#include <SImage.h>
#include <SImageIO.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <math.h>
#include <stdlib.h>


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
	return (-(pow(red_val - m[0], 2) / (pow(sig[0], 2)))
			* -(pow(green_val - m[1], 2) / (pow(sig[1], 2)))
			* -(pow(blue_val - m[2], 2) / (pow(sig[2], 2))));
}





void compute_unary_cost(vector<SDoublePlane> &D, SDoublePlane left_image, SDoublePlane right_image, int w, int d)
{
    double sum = 0;

    for(int i=w; i<left_image.rows()-w; i++)
    {
	for(int j=w; j<left_image.cols()-w; j++)
	{
	    for(int u=-w; u<w; u++)
	    {
		for(int v=-w; v<w; v++)
		{
		    sum += pow(left_image[i+u][j+v] - right_image[i+u][j+v+d], 2);
		}
	    }
	    D[d][i][j] = sum;
	    sum = 0;
	}
    }
}




SDoublePlane mrf_stereo(const SDoublePlane &left_image, const SDoublePlane &right_image)
{
  // implement this in step 4...
  //  this placeholder just returns a random disparity map
  SDoublePlane result(left_image.rows(), left_image.cols());

  for(int i=0; i<left_image.rows(); i++)
    for(int j=0; j<left_image.cols(); j++)
      result[i][j] = rand() % 256;

  return result;
}



int main(int argc, char *argv[])
{
  if(argc != 4 && argc != 3)
    {
      cerr << "usage: " << argv[0] << " image_file1 image_file2 [gt_file]" << endl;
      return 1;
    }

  string input_filename1 = argv[1], input_filename2 = argv[2];
  string gt_filename;
  if(argc == 4)
    gt_filename = argv[3];

  // read in images and gt
  SDoublePlane image1 = SImageIO::read_png_file(input_filename1.c_str());
  SDoublePlane image2 = SImageIO::read_png_file(input_filename2.c_str());
  SDoublePlane gt;
  if(gt_filename != "")
  {
    gt = SImageIO::read_png_file(gt_filename.c_str());
    // gt maps are scaled by a factor of 3, undo this...
    for(int i=0; i<gt.rows(); i++)
      for(int j=0; j<gt.cols(); j++)
        gt[i][j] = gt[i][j] / 3.0;
  }
  
  // do stereo using mrf
  SDoublePlane disp3 = mrf_stereo(image1, image2);
  SImageIO::write_png_file("disp_mrf.png", disp3, disp3, disp3);

  // Measure error with respect to ground truth, if we have it...
  if(gt_filename != "")
    {
      double err=0;
      for(int i=0; i<gt.rows(); i++)
	for(int j=0; j<gt.cols(); j++)
	  err += sqrt((disp3[i][j] - gt[i][j])*(disp3[i][j] - gt[i][j]));

      cout << "MRF stereo technique mean error = " << err/gt.rows()/gt.cols() << endl;

    }

  return 0;
}

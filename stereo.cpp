#include <SImage.h>
#include <SImageIO.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <math.h>

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

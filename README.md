MRF
===
In this assignment you’ll implement two different computer vision algorithms that can be posed using the
same underlying mathematical and algorithmic framework: Markov Random Fields (MRFs).
The first algorithm is for semi-supervised segmentation. The program will take a color image and a set of
“seed points” that indicate where some of the foreground and background pixels in the image are located.
These seed points could be quick strokes drawn manually by a human, like in Figure 1. The program will
then use these seeds in an MRF model to produce a complete segmentation of the image, partitioning the
image into foreground and background (Figure 1). The second algorithm is for resolving stereo: given two
images of the same scene taken from two different camera angles, infer the depth of each pixel by computing
disparities between the two images.
Despite the fact that these two algorithms seem quite different, they can both be posed as MRF inference
problems “under the hood,” so that much of the code can be shared between the two programs.

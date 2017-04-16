#include <opencv2/opencv.hpp>
#include<opencv2/core/core.hpp>
#include<algorithm>
#include "EST_El.h"
#include "Ellipse.h"
#include "Tangent.h"
#include "est.h"

int eight_neighbor(const cv::Point2i &ce, int x, int y);

void vote_theta(std::vector<cv::Point2f> &ce, cv::Point2i P, int *vote);
void vote_rad(std::vector<float> radi, std::vector<int> &r);
void detect_line(int th, int r);
void M_estimate(cv::Point2f &e, cv::Point2i c);
void colorExtraction(cv::Mat* src, cv::Mat* dst,
	int code,
	int ch1Lower, int ch1Upper, //@comment H(色相)　最小、最大
	int ch2Lower, int ch2Upper, //@comment S(彩度)　最小、最大
	int ch3Lower, int ch3Upper  //@comment V(明度)　最小、最大
	);
float calculate_distance(cv::Point2f, cv::Point2i);
float diff_Cg(float g, int x);
float diff_Cf(float f, int y);
float calc_Wf(cv::Point2f p0, cv::Point2i p1, float C);
float calc_Wg(cv::Point2f p0, cv::Point2i p1, float C);
void process_4(cv::Point2i O, std::vector<EST_El> &El);
std::vector<cv::Point3i> Kimura(cv::Mat* img, std::vector<cv::Point2i> &e, int *vote_k);
std::vector<cv::Point2i> vote_center(std::vector<std::vector<int>> &arr, std::vector<EST_El> &El);

std::vector<El> Hough_Ellipse(std::vector<cv::Point2f> HC, cv::Mat* img);
void Estimate_Ellipse(std::vector<EdgeInf> &EI, std::vector<EST_El> &Es);
void Hough_EllipseD(std::vector<cv::Point2i> &HC, std::vector<El> &e, cv::Mat* img);
bool cmp(cv::Point3i A, cv::Point3i B);
void RemoveEdge(cv::Mat *img);
void RemoveEdge2(cv::Mat *img);



void OnImageProcessing(cv::Mat3b* src, cv::Mat1b*, cv::Mat1b*);
void storeEdges(std::vector<cv::Point2i> &, cv::Mat1b *);
std::vector<Tangent>  deriveTangent(std::vector<Tangent> &Tan, std::vector<cv::Point2i> &, cv::Mat1b *);
void derive_Theta_R(std::vector<cv::Point2i> &, cv::Point2i &, int *, std::vector<int> &);
void calculate_center(cv::Point2i *c, cv::Point2i *e1, cv::Point2i *e2, int *vote, std::vector<int> &rad);

double get_vector_length(cv::Point2i v);
double dot_product(cv::Point2i vl, cv::Point2i vr);
double AngleOf2Vector(cv::Point2i A, cv::Point2i B);

void calculate_theta(cv::Point2f *c, cv::Point2f *e1, cv::Point2f *e2, int *vote);
void est_EllipseCenter(std::vector<Tangent> &Tan, cv::Mat1b *bin_img, std::vector<cv::Point2i> &C);




cv::Rect calcArea(std::vector<cv::Point2i> &points);
inline bool cmp2ix(cv::Point2i A, cv::Point2i B);
inline bool cmp2iy(cv::Point2i A, cv::Point2i B);


void vote_Center(std::vector<EST_El> &El, cv::Mat1b *bin, std::vector<cv::Point2i> &C);

void dst_center(std::vector<cv::Point2i> &c_set, std::vector<EST_El> &El, std::vector<cv::Point2i> &C);


void RemoveEdge3(cv::Mat *img);

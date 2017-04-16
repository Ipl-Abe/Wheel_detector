#include<opencv2/core/core.hpp>

class Tangent{
private:
	cv::Point2i contact;
	int th;
	int r;
public:
	Tangent(cv::Point2i contact, int th, int r);
	int get_th();
	int get_r();
	cv::Point2i get_contact();
};

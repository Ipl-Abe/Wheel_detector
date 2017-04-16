#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>
class ES{
private:
	cv::Point2f Pi, Pj, Qij, M;
	float a;
	float b, Pi_th, Pj_th;
public:
	ES(cv::Point2f Pi, cv::Point2f Pj, cv::Point2f Qij);
	void setPi_th(float Pi_th);
	void setPj_th(float Pj_th);
	float getDPi();
	float getDPj();
	cv::Point2f get_Pi();
	cv::Point2f get_Pj();
	cv::Point2f get_Qij();
	cv::Point2f get_M();
	void calculate_M(cv::Point2f Pi, cv::Point2f Pj);
	void set_l(cv::Point2f Q, cv::Point2f M);
	float get_a();
	float get_b();
};
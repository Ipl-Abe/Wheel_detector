#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>



class EST_El{
private:
	cv::Point2i Pi, Pj, Qij, M;
	float a;
	int b, Pi_th, Pj_th;
public:
	EST_El(cv::Point2i Pi, cv::Point2i Pj, cv::Point2i Qij);
	void setPi_th(int Pi_th);
	void setPj_th(int Pj_th);
	float getDPi();
	float getDPj();
	cv::Point2i get_Pi();
	cv::Point2i get_Pj();
	cv::Point2i get_Qij();
	cv::Point2i get_M();
	void calculate_M(cv::Point2i Pi, cv::Point2i Pj);
	void set_l(cv::Point2i Q, cv::Point2i M);
	float get_a();
	int get_b();
};








class EdgeInf{
private:
	cv::Point2i E, C;
	int theta;
	int ID;
public:
	EdgeInf(int ID, cv::Point2i E, cv::Point2i C, int theta);
	int get_ID();
	cv::Point2i get_E();
	cv::Point2i get_C();
	int get_theta();
};

class List{
private:
	cv::Point2i O;

public:
	List(cv::Point2i O);

};
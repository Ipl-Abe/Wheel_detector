#include "est.h"

ES::ES(cv::Point2f Pi, cv::Point2f Pj, cv::Point2f Qij){
	this->Pi = Pi;
	this->Pj = Pj;
	this->Qij = Qij;
}
void ES::setPi_th(float Pi_th){
	this->Pi_th = Pi_th;
}
void ES::setPj_th(float Pj_th){
	this->Pj_th = Pj_th;
}
float  ES::getDPi(){
	return -(cos(this->Pi_th*CV_PI / 180) / sin(this->Pi_th*CV_PI / 180));
}

float  ES::getDPj(){

	return  -(cos(this->Pj_th*CV_PI / 180) / sin(this->Pj_th*CV_PI / 180));
}
cv::Point2f ES::get_Pi(){
	return Pi;
}

cv::Point2f ES::get_Pj(){
	return Pj;
}

cv::Point2f ES::get_Qij(){
	return Qij;
}
cv::Point2f ES::get_M(){
	return M;
}

void ES::calculate_M(cv::Point2f Pi, cv::Point2f Pj){
	this->M.x = (Pi.x + Pj.x) / 2;
	this->M.y = (Pi.y + Pj.y) / 2;
}
void ES::set_l(cv::Point2f Q, cv::Point2f M){
	if (M.x - Q.x == 0.0) this->a = 0;
	else this->a = (float)(M.y - Q.y) / (M.x - Q.x);
	this->b = M.y - a*M.x;
}
float ES::get_a(){
	return a;
}
float ES::get_b(){
	return b;
}



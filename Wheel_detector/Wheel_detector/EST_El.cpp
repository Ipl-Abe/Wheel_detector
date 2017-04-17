#include "EST_El.h"



EST_El::EST_El(cv::Point2i Pi, cv::Point2i Pj, cv::Point2i Qij){
	this->Pi = Pi;
	this->Pj = Pj;
	this->Qij = Qij;
}
void EST_El::setPi_th(int Pi_th){
	this->Pi_th = Pi_th;
}
void EST_El::setPj_th(int Pj_th){
	this->Pj_th = Pj_th;
}
float  EST_El::getDPi(){
	return -(cos(this->Pi_th*CV_PI / 180) / sin(this->Pi_th*CV_PI / 180));
}

float  EST_El::getDPj(){

	return  -(cos(this->Pj_th*CV_PI / 180) / sin(this->Pj_th*CV_PI / 180));
}
cv::Point2i EST_El::get_Pi(){
	return Pi;
}

cv::Point2i EST_El::get_Pj(){
	return Pj;
}

cv::Point2i EST_El::get_Qij(){
	return Qij;
}
cv::Point2i EST_El::get_M(){
	return M;
}

void EST_El::calculate_M(cv::Point2i Pi, cv::Point2i Pj){
	this->M.x = (Pi.x + Pj.x) / 2;
	this->M.y = (Pi.y + Pj.y) / 2;
}
void EST_El::set_l(cv::Point2i Q, cv::Point2i M){
	if (M.x - Q.x == 0.0) this->a = 0;
	else this->a = (float)(M.y - Q.y) / (M.x - Q.x);
	this->b = M.y - a*M.x;
}
float EST_El::get_a(){
	return a;
}
int EST_El::get_b(){
	return b;
}




EdgeInf::EdgeInf(int ID, cv::Point2i E, cv::Point2i C, int theta){
	this->ID = ID;
	this->E = E;
	this->C = C;
	this->theta = theta;
}

int EdgeInf::get_ID(){
	return ID;
}

cv::Point2i EdgeInf::get_E(){
	return E;
}
cv::Point2i EdgeInf::get_C(){
	return C;
}

int EdgeInf::get_theta(){
	return theta;
}


List::List(cv::Point2i O){
	this->O = O;
}
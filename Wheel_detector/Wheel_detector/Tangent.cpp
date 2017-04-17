#include "Tangent.h"

Tangent::Tangent(cv::Point2i contact, int th, int r){
	this->contact = contact;
	this->th = th;
	this->r = r;
}

int Tangent::get_th(){
	return th;
}

int Tangent::get_r(){
	return r;
}

cv::Point2i Tangent::get_contact(){
	return contact;
}
#include "Ellipse.h"

El::El(int x, int y, int th, int a, int b){
	this->x = x;
	this->y = y;
	this->th = th;
	this->a = a;
	this->b = b;

}

int El::get_x(){
	return x;
}
int El::get_y(){
	return y;
}
int El::get_th(){
	return th;
}
int El::get_a(){
	return a;
}
int El::get_b(){
	return b;
}
int El::get_count(){
	return count;
}
int El::set_count(int num){
	count = num;
	return 1;
}


El_para::El_para(cv::Point2i O){

}
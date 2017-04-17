#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/opencv.hpp>
#include <iostream>
#include<vector>
#include <math.h>
#include <algorithm>
#include <stdlib.h>
#include "ImgProcess.h"




void OnImage(){
cv::Mat1b Mask; //�}�X�N�摜�p�ϐ�
std::string img_Name = "../snap3.png"; //���͉摜�p���O
cv::Mat3b img = cv::imread(img_Name);  //�摜�ǂݍ���
cv::Mat1b stes; 
cv::Mat1b canny;  //�L���j�[�摜�p�ϐ�

if (!img.data){
	exit(0);
}

cv::resize(img,img,cv::Size(),0.05,0.05); //�摜�̃��T�C�Y

colorExtraction(&img, &Mask, CV_BGR2HSV, 0, 180, 70, 255, 70, 255); //�A�[���F���o
cv::Mat1b gray;	//�O���[�X�P�[���p�ϐ�
cv::cvtColor(img, gray, CV_BGR2GRAY); //�O���[�摜�ϐ�

cv::GaussianBlur(gray, gray, cv::Size(3, 3), 2, 2); //�K�E�V�A���t�B���^

int d_value = 0;  //�F���o�p�����[�^1
int c_value = 0;  //�F���o�p�����[�^2

cv::namedWindow("colorExt", 1);
cv::createTrackbar("c1", "colorExt", &d_value, 255);
cv::createTrackbar("c2", "colorExt", &c_value, 255);

//while (1){
//	cv::Canny(gray, canny, d_value,c_value, 3, false);//@comment canny�t�B���^
//	imshow("test_canny",canny);
//	if (cv::waitKey(10) == 'q')
//	{
//		break;
//	}
//}


cv::Canny(gray, canny, 97, 166, 3, false);//@comment canny�t�B���^
canny.copyTo(stes, Mask);

cv::imshow("cany", canny);
cv::imshow("ca", stes);
cv::waitKey(100);

OnImageProcessing(&img, &stes, &canny);
cv::imshow("fin_sss", canny);
std::cout << "end of program" << std::endl;
cv::waitKey(0);


}

int main(){

	OnImage();
	return 0;
}
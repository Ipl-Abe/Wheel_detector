#include "ImgProcess.h"


cv::Mat element = cv::Mat::ones(3,3,CV_8UC1); //dilate�����ɕK�v�ȍs��


void OnImageProcessing(cv::Mat3b* src, cv::Mat1b* bin_img, cv::Mat1b* bin_origin){
			
}





void colorExtraction(cv::Mat* src, cv::Mat* dst,
	int code,
	int ch1Lower, int ch1Upper, //@comment H(�F��)�@�ŏ��A�ő�
	int ch2Lower, int ch2Upper, //@comment S(�ʓx)�@�ŏ��A�ő�
	int ch3Lower, int ch3Upper  //@comment V(���x)�@�ŏ��A�ő�
	)
{
	cv::Mat colorImage;
	int lower[3];
	int upper[3];

	cv::Mat lut = cv::Mat(256, 1, CV_8UC3);

	cv::cvtColor(*src, colorImage, code);

	lower[0] = ch1Lower;
	lower[1] = ch2Lower;
	lower[2] = ch3Lower;

	upper[0] = ch1Upper;
	upper[1] = ch2Upper;
	upper[2] = ch3Upper;

	for (int i = 0; i < 256; i++) {
		for (int k = 0; k < 3; k++) {
			if (lower[k] <= upper[k]) {
				if ((lower[k] <= i) && (i <= upper[k])) {
					lut.data[i*lut.step + k] = 255;
				}
				else {
					lut.data[i*lut.step + k] = 0;
				}
			}
			else {
				if ((i <= upper[k]) || (lower[k] <= i)) {
					lut.data[i*lut.step + k] = 255;
				}
				else {
					lut.data[i*lut.step + k] = 0;
				}
			}
		}
	}
	//@comment LUT���g�p���ē�l��
	cv::LUT(colorImage, lut, colorImage);

	//@comment Channel���ɕ���
	std::vector<cv::Mat> planes;
	cv::split(colorImage, planes);

	//@comment �}�X�N���쐬
	cv::Mat maskImage;
	cv::bitwise_and(planes[0], planes[1], maskImage);

	cv::bitwise_and(maskImage, planes[2], maskImage);
	//dilate(maskImage, maskImage, element, cv::Point(-1, -1), 3);
	dilate(maskImage, maskImage, element, cv::Point(-1, -1), 3);
	erode(maskImage, maskImage, element, cv::Point(-1, -1), 1);
	cv::imshow("test", ~maskImage);
	cv::waitKey(100);
	//@comemnt �o��
	cv::Mat maskedImage;
	//src->copyTo(maskedImage, ~maskImage);
	//*dst = maskedImage;
	*dst = ~maskImage;
}
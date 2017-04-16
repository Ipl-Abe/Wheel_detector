#include "ImgProcess.h"


cv::Mat element = cv::Mat::ones(3,3,CV_8UC1); //dilate処理に必要な行列


void OnImageProcessing(cv::Mat3b* src, cv::Mat1b* bin_img, cv::Mat1b* bin_origin){
	//cv::Mat s = cv::imread("../snap3.png", 1);
	//cv::Mat r;
	//s.copyTo(origin);
	//cv::resize(s, s, cv::Size(), 0.05, 0.05);
	//s.copyTo(test_img);

	cv::Mat1b Bin;


	/*************エッジ点の格納を行う**************/
	std::vector<cv::Point2i> EdgeArray;
	storeEdges(EdgeArray, bin_img); //エッジ点を格納
	std::cout <<"Total Edge number : "<< EdgeArray.size() << std::endl; //格納されたエッジ点数を表示



	/*****あるエッジ点を注目点としてその周り13x13のうちの2ペアEi,Ejを設定し、
										それらで注目エッジ点の接線の方程式を求める****************/
	std::vector<Tangent> TV;

	std::cout << "start deriveTangent " << std::endl;

	std::vector<Tangent> tanVec = deriveTangent(TV, EdgeArray, bin_img);

	std::cout << "end deriveTangent " << std::endl;



}





void colorExtraction(cv::Mat* src, cv::Mat* dst,
	int code,
	int ch1Lower, int ch1Upper, //@comment H(色相)　最小、最大
	int ch2Lower, int ch2Upper, //@comment S(彩度)　最小、最大
	int ch3Lower, int ch3Upper  //@comment V(明度)　最小、最大
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
	//@comment LUTを使用して二値化
	cv::LUT(colorImage, lut, colorImage);

	//@comment Channel毎に分解
	std::vector<cv::Mat> planes;
	cv::split(colorImage, planes);

	//@comment マスクを作成
	cv::Mat maskImage;
	cv::bitwise_and(planes[0], planes[1], maskImage);

	cv::bitwise_and(maskImage, planes[2], maskImage);
	//dilate(maskImage, maskImage, element, cv::Point(-1, -1), 3);
	dilate(maskImage, maskImage, element, cv::Point(-1, -1), 3);
	erode(maskImage, maskImage, element, cv::Point(-1, -1), 1);
	cv::imshow("test", ~maskImage);
	cv::waitKey(100);
	//@comemnt 出力
	cv::Mat maskedImage;
	//src->copyTo(maskedImage, ~maskImage);
	//*dst = maskedImage;
	*dst = ~maskImage;
}


void storeEdges(std::vector<cv::Point2i> &EdgeArray, cv::Mat1b * img){
	
	for (int j = 0; j < img->rows; j++){
		unsigned char *ptr = img->ptr<unsigned char>(j);
		for (int i = 0; i < img->cols; i++){
			std::cout << " " << (int)ptr[i] / 250;  //コンソール表示の際に
			if ((int)ptr[i] > 0 && j > 6 && j < img->rows - 6
				&& i >   6 && i < img->cols - 6){
				cv::Point2i edge_xy(i, j);
				EdgeArray.push_back(edge_xy);
			}
		}
		std::cout << std::endl;
	}
	

}


//接線の情報を取得
std::vector<Tangent> deriveTangent(std::vector<Tangent> &Tan, std::vector<cv::Point2i> &E, cv::Mat1b *img){
	std::vector<int> rad;
	std::vector<cv::Point2i>::iterator it;
	std::vector<cv::Point2i> EdgeSet;
	double s = 0;
	cv::Point2i AE, AB, e1, e2;
	cv::Point2i ate;

	for (auto it : E){
		int vote_k[361] = { 0 };
		int r_array[500] = { 0 };

		for (int j = it.y - 6; j < it.y + 7; j++){
			unsigned char *ptr = img->ptr<unsigned char>(j);
			for (int i = it.x - 6; i < it.x + 7; i++){
				//@comment 中心と同じ場合、画像の領域から外れる場合はpass
				if ((j != it.y && i != it.x) &&
					(j > 0 && j < img->rows) &&
					(i > 0 && i < img->cols)
					){
					if ((int)ptr[i] > 0){
						cv::Point2i ce;
						ce.x = it.x;
						ce.y = it.y;
						//@comment 8近傍に存在するエッジは省く。
						if (!eight_neighbor(ce, i, j)){
							cv::Point2i cc;
							cc.x = i;
							cc.y = j;
							EdgeSet.push_back(cc);
						}
					}
				}
			}
		}

		cv::Point2i attentionEdge(it.x, it.y);
		derive_Theta_R(EdgeSet, attentionEdge, vote_k, rad);


		int max_count = 1;
		int index = 0;
		int start = 0;

		for (int i = 0; i < 361; i++){
			if (max_count <= vote_k[i]){
				if (max_count == vote_k[i]){
					index = (start + i) / 2;
				}
				else{
					max_count = vote_k[i];
					index = i;
					start = i;
				}
			}
		}
		std::vector<int>::iterator itr;
		int max_rcount = 1;
		int r_start = 0;
		size_t r_index = 0;
		int r_count = 0;
		for (auto itr : rad){
			if (max_rcount <= itr)
			{
				if (max_rcount == itr){
					r_index = (r_start + r_count) / 2;
				}
				else{
					max_rcount = itr;
					r_index = r_count;
					r_start = r_count;
				}
			}
			r_count++;
		}

		if (index != 0 && r_index != 0){
			Tangent t = Tangent(attentionEdge, index - 180, r_index);
			Tan.push_back(t);
		}
		EdgeSet.clear();
	}
	return Tan;
}
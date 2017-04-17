#include "ImgProcess.h"


cv::Mat element = cv::Mat::ones(3,3,CV_8UC1); //dilate処理に必要な行列

int vote_k[181] = { 0 };
int r_array[500] = { 0 };

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


	std::vector<cv::Point2i> C;
	est_EllipseCenter(src,tanVec, bin_img, C);



	/////////////////////////////////////////////////////////
	/////重複要素削除
	/////////////////////////////////////////////////////////
	std::sort(C.begin(), C.end(), cmp2iy);
	C.erase(std::unique(C.begin(), C.end()), C.end());
	std::sort(C.begin(), C.end(), cmp2ix);
	C.erase(std::unique(C.begin(), C.end()), C.end());
	std::sort(C.begin(), C.end(), cmp2iy);
	C.erase(std::unique(C.begin(), C.end()), C.end());
	std::sort(C.begin(), C.end(), cmp2ix);
	C.erase(std::unique(C.begin(), C.end()), C.end());
	////////////////////////////////////////////////////
	////////////////////////////////////////////////////

	std::vector<El> v;

	std::cout << "get center number : " << C.size() << std::endl;
	//Hough_EllipseD(C,v,bin_img);
	Hough_EllipseD(src,C, v, bin_origin);


	std::cout << "fin the program" << std::endl;

	cvWaitKey(0);


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


void storeEdges(std::vector<cv::Point2i> &EdgeArray, cv::Mat1b * bin_img){
	
	for (int j = 0; j < bin_img->rows; j++){
		unsigned char *ptr = bin_img->ptr<unsigned char>(j);
		for (int i = 0; i < bin_img->cols; i++){
			std::cout << " " << (int)ptr[i] / 250;  //コンソール表示の際に
			//ピクセル値が0より大きい場合 && 画像の6ピクセル内側に内側にピクセルが
			//存在していればそのピクセルを格納する
			if ((int)ptr[i] > 0 &&    
				6 < j && j < bin_img->rows - 6 &&
				6 < i && i < bin_img->cols - 6){
				cv::Point2i edgeCoordinate_xy(i, j);
				EdgeArray.push_back(edgeCoordinate_xy);
			}
		}
		std::cout << std::endl;//表示のための改行
	}
	

}


//接線の情報を取得
std::vector<Tangent> deriveTangent(std::vector<Tangent> &Tan, std::vector<cv::Point2i> &E, cv::Mat1b *img){
	std::vector<int> rad; //円の半径格納用ベクタ配列
	std::vector<cv::Point2i>::iterator it;
	std::vector<cv::Point2i> EdgeSet;
	double s = 0;
	cv::Point2i AE, AB, e1, e2;
	cv::Point2i ate;


	//格納されたすべてのエッジに対して処理を行う
	for (auto it : E){
		int vote_k[361] = { 0 };
		int r_array[500] = { 0 };

		//注目エッジ点に対して13x13の領域を見る
		for (int j = it.y - 6; j < it.y + 7; j++){
			unsigned char *ptr = img->ptr<unsigned char>(j);
			for (int i = it.x - 6; i < it.x + 7; i++){
				//@comment 中心と同じ場合、画像の領域から外れる場合はpass
				if ((j != it.y && i != it.x) &&
					(j > 0 && j < img->rows) &&
					(i > 0 && i < img->cols)
					){
					if ((int)ptr[i] > 0){
						cv::Point2i near_edge = cv::Point2i(it.x, it.y);
						//ce.x = it.x;
						//ce.y = it.y;
						//@comment 8近傍に存在するエッジは省く。
						if (!eight_neighbor(near_edge, i, j)){
							cv::Point2i cc = cv::Point2i(i,j);
							//cc.x = i;
							//cc.y = j;
							EdgeSet.push_back(cc);
						}
					}
				}
			}
		}

		cv::Point2i attentionEdge(it.x, it.y);
		derive_Theta_R(EdgeSet, attentionEdge, vote_k, rad); //円の半径と角度を投票によって求める


		int max_count = 1;
		int index = 0;
		int start = 0;

		//もっともカウントが大きい角度をその円の角度とする
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
		//最もカウントが大きい半径をその円の半径とする
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
		//Tangentクラスに注目エッジの情報を登録する
		if (index != 0 && r_index != 0){
			Tangent t = Tangent(attentionEdge, index - 180, r_index);
			Tan.push_back(t);
		}
		EdgeSet.clear();
	}
	return Tan;
}


//@comment 8近傍の存在をチェック
inline int eight_neighbor(const cv::Point2i &ce, int x, int y){
	if (ce.x == x && ce.y == y - 1 ||
		ce.x == x && ce.y == y + 1 ||
		ce.x == x - 1 && ce.y == y - 1 ||
		ce.x == x - 1 && ce.y == y + 1 ||
		ce.x == x + 1 && ce.y == y - 1 ||
		ce.x == x + 1 && ce.y == y + 1 ||
		ce.x == x + 1 && ce.y == y ||
		ce.x == x - 1 && ce.y == y
		){
		return 1;
	}
	return 0;
}
//ベクトルの長さを計算する
inline double get_vector_length(cv::Point2i v) {
	return pow((v.x * v.x) + (v.y * v.y), 0.5);
}

//ベクトル内積
inline double dot_product(cv::Point2i vl, cv::Point2i vr) {
	return vl.x * vr.x + vl.y * vr.y;
}
//２つのベクトルABのなす角度θを求める
inline double AngleOf2Vector(cv::Point2i A, cv::Point2i B)
{
	//　※ベクトルの長さが0だと答えが出ませんので注意してください。

	//ベクトルAとBの長さを計算する
	double length_A = get_vector_length(A);
	double length_B = get_vector_length(B);

	//内積とベクトル長さを使ってcosθを求める
	double cos_sita = dot_product(A, B) / (length_A * length_B);

	//cosθからθを求める
	double sita = acos(cos_sita);

	//ラジアンでなく0～180の角度でほしい場合はコメント外す
	sita = sita * 180.0 / CV_PI;

	return sita;
}

//あるエッジ点でのθとRを投票によって求める
void derive_Theta_R(std::vector<cv::Point2i> &edges, cv::Point2i &ate, int *vote_k, std::vector<int> &rad){
	std::vector<cv::Point2i>::iterator vb, ve; //配列用イテレータ
	cv::Point2i e1, e2, center = ate;
	cv::Point2i AE, AB;
	double s = 0;

	//@comment すべてのエッジ組み合わせに対して処理を行う
	for (auto vb : edges){
		for (auto ve : edges){
			//@comment 同じエッジ点は処理しない
			if (!(vb.x == ve.x && vb.y == ve.y)){
				e1.x = vb.x;
				e1.y = vb.y;
				e2.x = ve.x;
				e2.y = ve.y;
				AE.x = (e2.x - ate.x);
				AE.y = (e2.y - ate.y);
				AB.x = (e1.x - ate.x);
				AB.y = (e1.y - ate.y);
				s = AngleOf2Vector(AB, AE);
				//ate(attention_edge)とPi,Pjのなす角が６０度以下の場合は処理を行わない。
				if (s >= 60){
					calculate_center(&center, &e1, &e2, vote_k, rad);
				}
			}
		}
	}
}



//@comment 与えられた３点から円の中心を推定する
// ピクセルの量子誤差を考慮してピクセルの４隅を設定する
void calculate_center(cv::Point2i *c, cv::Point2i *e1, cv::Point2i *e2, int *vote, std::vector<int> &rad){
	float r = 0;
	std::vector<cv::Point2f> ej, ek;
	std::vector<float> radius;
	cv::Point2f e11, e12, e13, e14;
	cv::Point2f e21, e22, e23, e24;
	std::vector<cv::Point2f>::iterator it1, it2;
	cv::Mat A, B, C;
	float a11 = 0;
	float a12 = 0;
	float a21 = 0;
	float a22 = 0;
	float b11 = 0;
	float b21 = 0;
	float x = 0;
	float y = 0;
	std::vector<cv::Point2f> center;
	cv::Mat_<float>::iterator it_x, it_y;
	cv::Point2f cpara;


	e11.x = e1->x - 0.5;
	e11.y = e1->y - 0.5;
	e12.x = e1->x + 0.5;
	e12.y = e1->y - 0.5;
	e13.x = e1->x + 0.5;
	e13.y = e1->y + 0.5;
	e14.x = e1->x - 0.5;
	e14.y = e1->y + 0.5;

	e21.x = e2->x - 0.5;
	e21.y = e2->y - 0.5;
	e22.x = e2->x + 0.5;
	e22.y = e2->y - 0.5;
	e23.x = e2->x + 0.5;
	e23.y = e2->y + 0.5;
	e24.x = e2->x - 0.5;
	e24.y = e2->y + 0.5;

	ej.push_back(e11);
	ej.push_back(e12);
	ej.push_back(e13);
	ej.push_back(e14);
	ek.push_back(e21);
	ek.push_back(e22);
	ek.push_back(e23);
	ek.push_back(e24);

	//3点を通る円の中心を求める
	for (auto it1 : ej){
		for (auto it2 : ek){
			a11 = 2 * (c->x - it1.x);
			a12 = 2 * (c->y - it1.y);
			a21 = 2 * (it2.x - c->x);
			a22 = 2 * (it2.y - c->y);
			b11 = (c->x * c->x) + (c->y*c->y) - (it1.x*it1.x) - (it1.y * it1.y);
			b21 = (it2.x * it2.x) + (it2.y*it2.y) - (c->x*c->x) - (c->y * c->y);
			A = (cv::Mat_<float>(2, 2) << a11, a12, a21, a22);
			B = (cv::Mat_<float>(2, 1) << b11, b21);
			C = (cv::Mat_<float>(2, 1) << x, y);
			C = A.inv() * B;//行列演算により座標取得

			it_x = C.begin<float>();
			it_y = C.end<float>() - 1;
			cpara.x = *it_x;
			cpara.y = *it_y;
			center.push_back(cpara);
			//円の半径取得
			r = sqrt(((c->x - (*it_x)) *(c->x - (*it_x))) + ((c->y - (*it_y))*(c->y - (*it_y))));
			radius.push_back(r);
		}
	}


	vote_theta(center, *c, vote);//角度に関する投票

	vote_rad(radius, rad);// 半径に対する投票
}


void vote_theta(std::vector<cv::Point2f> &ce, cv::Point2i P, int *vote){
	std::vector<cv::Point2f>::iterator it;
	std::vector<float> th;
	int max_theta = 0;
	int min_theta = 180;
	int max_theta2 = -180;
	int min_theta2 = 0;
	float theta = 0;

	for (auto it : ce){
		if ((it.x - P.x)> 0 && (it.y - P.y)< 0){ //第一象限
			theta = -atan2(it.y - P.y, it.x - P.x) * 180 / CV_PI;
		}
		else if ((it.x - P.x)< 0 && (it.y - P.y) < 0){//第二象限
			theta = -atan2(it.y - P.y, it.x - P.x) * 180 / CV_PI;
			if (theta < 0){
			}
		}
		else if ((it.x - P.x) > 0 && (it.y - P.y) > 0){//第四象限
			theta = -atan2(it.y - P.y, it.x - P.x) * 180 / CV_PI;
		}
		else if ((it.x - P.x) < 0 && (it.y - P.y) > 0){//第三象限
			theta = -atan2(it.y - P.y, it.x - P.x) * 180 / CV_PI;
		}
		th.push_back(theta);
	}
	std::vector<float>::iterator th_it;
	std::sort(th.begin(), th.end());
	for (auto th_it : th){
		if (th_it > 0){
			if (th_it > max_theta){
				max_theta = th_it;
			}
			if (th_it < min_theta){
				min_theta = th_it;
			}
		}
		else if (th_it < 0){
			if (th_it > max_theta2){
				max_theta2 = th_it;
			}
			if (th_it < min_theta2){
				min_theta2 = th_it;
			}
		}
	}
	if (max_theta - (max_theta2 + 10) < 0){
		if (max_theta != 0){
			min_theta = min_theta2;
		}
	}
	//vote_theta内
	for (int i = min_theta; i < max_theta; i++){
		*(vote + 180 + i) += 1;
	}
	for (int i = min_theta2; i < max_theta2; i++){
		*(vote + 180 + i) += 1;
	}
}
void vote_rad(std::vector<float> radi, std::vector<int> &r){
	//std::vector<float>::iterator so;
	std::sort(radi.begin(), radi.end()); //rを降順にソート
	std::vector<float>::iterator so1, so2;
	so1 = radi.begin();
	so2 = radi.end() - 1;
	float min_r = *so1;
	float max_r = *so2;

	std::vector<int>::iterator rr;
	rr = r.begin();
	for (int i = 0; i < r.size(); i++){
		r_array[i] = *rr;
		rr++;
	}
	std::vector<int> v1;
	for (int i = 0; i < max_r; i++){
		int num = 0;
		if (i >= min_r && i <= max_r){
			num = r_array[i] + 1;
		}
		v1.push_back(num);
	}
	r.clear();
	copy(v1.begin(), v1.end(), back_inserter(r));
}


inline bool cmp(cv::Point3i A, cv::Point3i B){

	return A.z > B.z ? 1 : 0;
}
inline bool cmp2ix(cv::Point2i A, cv::Point2i B){

	return A.x > B.x ? 1 : 0;
}
inline bool cmp2iy(cv::Point2i A, cv::Point2i B){

	return A.y > B.y ? 1 : 0;
}
void est_EllipseCenter(cv::Mat3b *src, std::vector<Tangent> &Tan, cv::Mat1b *bin_img, std::vector<cv::Point2i> &C){
	const int count = 0;
	cv::Mat te, tte;
	std::vector<Tangent>::iterator tan1, tan2;
	cv::Point2i e1, e2;
	int flag = 0;
	int i = 0;
	int y = 0;
	float co1 = 0;
	float co2 = 0;
	float si1 = 0;
	float si2 = 0;
	float R1 = 0;
	float R2 = 0;
	float a = 0;
	float b = 0;
	float c = 0;
	float d = 0;

	cv::Point2i Q, Q2, Q3;


	for (auto tan1 : Tan){
		int count = 0;
		for (auto tan2 : Tan){
			if (tan1.get_contact() != tan2.get_contact()){
				e1 = tan1.get_contact();
				e2 = tan2.get_contact();
				co1 = cos(tan1.get_th()*CV_PI / 180);
				co2 = cos(tan2.get_th()*CV_PI / 180);
				si1 = -sin(tan1.get_th()*CV_PI / 180);
				si2 = -sin(tan2.get_th()*CV_PI / 180);
				R1 = e1.x * co1 + e1.y *si1;
				R2 = e2.x * co2 + e2.y *si2;
				a = -(co1 / si1);
				b = R1 / si1;
				c = -(co2 / si2);
				d = R2 / si2;

				src->copyTo(te);
				//逆ハフ変換を用いた方法　接線
				for (i = e1.x - 20; i <e1.x + 50; i++){
					y = -(co1 / si1)* i + R1 / si1;
					//cout << "y: " << y << " x: "<< i << endl;
					circle(te, cv::Point(i, y), 1, cv::Scalar(0, 0, 255), -1, CV_AA);
				}
				for (i = e2.x - 20; i <e2.x + 50; i++){
					y = -(co2 / si2)* i + R2 / si2;
					//cout << "y: " << y << " x: "<< i << endl;
					cv::circle(te, cv::Point(i, y), 1, cv::Scalar(0, 255, 0), -1, CV_AA);
				}
				//cv::imshow("rete",te);
				//cv::waitKey(10);


				float dist = abs(sqrt((e1.x - e2.x) * (e1.x - e2.x) + (e1.y - e2.y)*(e1.y - e2.y)));

				//Qを求める

				Q.x = (b - d) / (c - a);
				Q.y = -(a*d - b*c) / (c - a);
				cv::circle(te, cv::Point(Q.x, Q.y), 3, cv::Scalar(255, 0, 0), -1, CV_AA);

				EST_El es = EST_El(e1, e2, Q);
				es.setPi_th(tan1.get_th());
				es.setPj_th(tan2.get_th());
				cv::Point2i M;
				es.calculate_M(es.get_Pi(), es.get_Pj());
				M = es.get_M();
				es.set_l(Q, M);
				cv::circle(te, cv::Point(M.x, M.y), 3, cv::Scalar(255, 0, 255), -1, CV_AA);
				//cv::imshow("circle", te);
				//cv::waitKey(10);
				float dist2 = abs(sqrt((M.x - Q.x) * (M.x - Q.x) + (M.y - Q.y)*(M.y - Q.y)));
				float dist3 = abs(e1.y - e2.y);
				float dist4 = abs(e1.x - e2.x);
				if ((dist > 10 &&
					dist2 > 4 &&
					dist3 > 2 &&
					dist4 > 2)
					//(Q.x > 0 && Q.x < te.rows) &&
					//(Q.y > 0 && Q.y < te.cols)
					){
					//e1.e2がある程度離れているならばそれらの間の３点目の接線を求める

					cv::Point2i e3;
					//int Y=0;
					float ta;
					ta = es.get_a();
					int B = es.get_b();

					for (i = M.x; i <Q.x; i++){
						y = ta * i + B;
						cv::circle(te, cv::Point(i, y), 1, cv::Scalar(0, 0, 0), -1, CV_AA);
					}
					//cv::imshow("ppap", te);
					//cv::waitKey(1);


					for (int n = M.x; n <= Q.x; n++){
						y = es.get_a() * n + es.get_b();

						if (y < 1 || y > bin_img->rows - 1 || n < 1 || n > bin_img->cols - 1){
							break;
						}
						//else if (((int)bin_img->at<unsigned char>(Y-1, n) > 0 ||
						//	(int)bin_img->at<unsigned char>(Y + 1, n) > 0 ||
						//	(int)bin_img->at<unsigned char>(Y, n-1) > 0 ||
						//	(int)bin_img->at<unsigned char>(Y, n+1) > 0 ||
						//	(int)bin_img->at<unsigned char>(Y-1, n-1) > 0 ||
						//	(int)bin_img->at<unsigned char>(Y - 1, n+1) > 0 ||
						//	(int)bin_img->at<unsigned char>(Y + 1, n - 1) > 0 ||
						//	(int)bin_img->at<unsigned char>(Y + 1, n + 1) > 0)){
						else if ((int)bin_img->at<unsigned char>(y, n) > 0){
							e3.x = n;
							e3.y = y;
							flag = 1;
							break;
						}
					}



					if (flag == 1){
						flag = 0;
						float e = 0;
						float f = 0;

						double ta_e1e2 = atan((double)(e1.y - e2.y) / (e1.x - e2.x)) * 180 / CV_PI;
						double co3 = cos(ta_e1e2 * CV_PI / 180);
						double si3 = -sin(ta_e1e2 * CV_PI / 180);

						float R3 = e3.x * co3 + e3.y *si3;
						e = -(co3 / si3);
						f = R3 / si3;
						src->copyTo(tte);
						for (int i = e3.x - 20; i <e3.x + 50; i++){
							y = -(co3 / si3)* i + R3 / si3;
							//cout << "y: " << y << " x: "<< i << endl;
							cv::circle(te, cv::Point(i, y), 1, cv::Scalar(200, 30, 255), -1, CV_AA);
						}

						//cv::imshow("popo",te);
						//cv::waitKey(1);
						Q2.x = (b - f) / (e - a);
						Q2.y = -(a*f - b*e) / (e - a);
						Q3.x = (d - f) / (e - c);
						Q3.y = -(c*f - d*e) / (e - c);

						//				//std::cout << " Q2 : " << Q2 << " Q3 : " << Q3 << std::endl;
						cv::circle(tte, cv::Point(Q2.x, Q2.y), 4, cv::Scalar(0, 0, 0), -1, CV_AA);
						cv::circle(tte, cv::Point(Q3.x, Q3.y), 4, cv::Scalar(10, 10, 10), -1, CV_AA);
						//cv::imshow("test",tte);
						//cv::waitKey(1);
						//
						//
						std::vector<cv::Point2i> points;
						points.push_back(Q);
						points.push_back(Q2);
						points.push_back(Q3);
						cv::Rect rcc = calcArea(points);
						//std::cout << rcc.x << ":" << rcc.width << ":" << rcc.y << ":" << rcc.height << std::endl;

						if ((Q2.x >= rcc.x && Q2.x <= rcc.width) &&
							(Q2.y >= rcc.y && Q2.y <= rcc.height) &&
							(Q3.x >= rcc.x && Q3.x <= rcc.width) &&
							(Q3.y >= rcc.y && Q3.x <= rcc.height)
							){

							std::vector<EST_El> VecEl;
							EST_El el1 = EST_El(e1, e2, Q);
							EST_El el2 = EST_El(e1, e3, Q2);
							EST_El el3 = EST_El(e2, e3, Q3);


							el1.setPi_th(tan1.get_th());
							el1.setPj_th(tan2.get_th());
							el1.calculate_M(e1, e2);
							el1.set_l(Q, el1.get_M());

							VecEl.push_back(el1);

							el2.setPi_th(tan1.get_th());
							el2.setPj_th(ta_e1e2);
							el2.calculate_M(e1, e3);
							el2.set_l(Q2, el2.get_M());

							VecEl.push_back(el2);

							el3.setPi_th(tan2.get_th());
							el3.setPj_th(ta_e1e2);
							el3.calculate_M(e2, e3);
							el3.set_l(Q3, el3.get_M());

							VecEl.push_back(el3);
							//
							//
							vote_Center(src, VecEl, bin_img, C);
							//
						}
						//				}
					}
				}
			}
		}
	}
}

cv::Rect calcArea(std::vector<cv::Point2i> &points){
	cv::Rect rc;
	int xmin = 0;
	int xmax = 0;
	int ymin = 0;
	int ymax = 0;
	std::vector<cv::Point2i>::iterator it;

	std::sort(points.begin(), points.end(), cmp2ix);
	it = points.begin();
	xmax = it->x;
	it = points.end() - 1;
	xmin = it->x;

	std::sort(points.begin(), points.end(), cmp2iy);
	it = points.begin();
	ymax = it->y;
	it = points.end() - 1;
	ymin = it->y;

	return rc = cv::Rect(xmin, ymin, xmax, ymax);
}

void vote_Center(cv::Mat3b* src, std::vector<EST_El> &El, cv::Mat1b *bin, std::vector<cv::Point2i> &C){
	cv::Mat cop;
	int varray[500][500] = { 0 };
	std::vector<EST_El>::iterator it;
	cv::Point2i Q;
	cv::Point2i M;
	int y = 0;
	int max_count = 0;
	cv::Point2i max_index;
	int tcount = 0;


	cv::Point2i ccc;
	std::vector<cv::Point2i> c_set;


	//std::cout << "in vote_center" << std::endl;
	for (auto it : El){
		src->copyTo(cop);
		Q = it.get_Qij();
		M = it.get_M();

		for (int i = 0; i <bin->cols; i++){
			y = it.get_a()*i + it.get_b();
			if (y < 0 || y > bin->rows) continue;
			varray[i][y]++;
			//cv::circle(cop, cv::Point(i, y), 1, cv::Scalar(200, 30, 255), -1, CV_AA);

		}
		//cv::imshow("vnj", cop);
		//cv::waitKey(20);

	}
	//for (int j = 0; j < bin->rows; j++){
	//	for (int i = 0; i < bin->cols; i++){
	//		std::cout << " " << varray[i][j];
	//		}
	//		std::cout << std::endl;
	//	}

	for (int j = 0; j < bin->rows; j++){
		for (int i = 0; i < bin->cols; i++){
			if (varray[i][j] >= 2){
				cv::Point2i ccc = cv::Point2i(i, j);
				ccc.x = i;
				ccc.y = j;
				c_set.push_back(ccc);
			}
		}
	}
	if (c_set.size() >= 1){
		dst_center(c_set, El, C);
	}
}



void  dst_center(std::vector<cv::Point2i> &c_set, std::vector<EST_El> &El, std::vector<cv::Point2i> &C){
	std::vector<cv::Point2i>::iterator cc;
	std::vector<EST_El>::iterator it;
	std::vector<cv::Point3i> vote;
	std::vector<cv::Point3i>::iterator vote_it;
	cv::Point3i ss;
	float a = 0;
	float b = 0;
	float d = 0;
	float res = 0;
	std::vector<float> bb;
	cv::Point2i center;

	for (auto cc : c_set){
		res = 0;
		for (auto it : El){
			a = it.get_a();

			d = abs(cc.y - cc.x * it.get_a() - it.get_b()) / sqrt(it.get_a()*it.get_a() + 1);
			res += d;
		}
		ss.x = cc.x;
		ss.y = cc.y;
		ss.z = res;
		vote.push_back(ss);

	}
	std::sort(vote.begin(), vote.end(), cmp);
	vote_it = vote.end() - 1;
	if (vote_it->z <= 0){
		std::cout << vote_it->z << std::endl;
		center.x = vote_it->x;
		center.y = vote_it->y;
		C.push_back(center);
	}
}



void Hough_EllipseD(cv::Mat3b* src, std::vector<cv::Point2i> &HC, std::vector<El> &v, cv::Mat *img){
	cv::Mat image2;

	int num = 0;
	int maxnum = 90;
	int mmx = 0;
	int px = 0;
	int py = 0;
	int Fx = 0, Fy = 0, Fx_ = 0, Fy_ = 0;
	int F = 0, P = 0;
	int X = 0, Y = 0;
	int rex = 0, rey = 0;
	int dist = 0;

	int count = 0;
	std::vector<cv::Point2i>::iterator it;

	for (it = HC.begin(); it != HC.end(); it++){
		mmx = 0;
		std::cout << "count" << count++ << std::endl;
		for (int l = 0; l < 180; l += 2){
			for (int m = src->cols / 8; m < src->cols / 3; m += 1){
				for (int n = src->rows / 8; n < src->rows / 2; n += 1){
					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
					El el = El::El(it->x, it->y, l, m, n);
					for (int o = 0; o < 360; o += 1){
						X = (int)(m)*cos(o*CV_PI / 180);
						Y = (int)(n)*sin(o*CV_PI / 180);
						rex = (int)X*cos(l*CV_PI / 180) - Y*sin(l*CV_PI / 180) + it->x;
						rey = (int)X*sin(l*CV_PI / 180) + Y*cos(l*CV_PI / 180) + it->y;
						if (rex >= 1 && rex < src->cols - 1 && rey >= 1 && rey < src->rows - 1){
							if ((int)img->at<unsigned char>(rey, rex) > 0){
								num++;
							}
						}
					}
					src->copyTo(image2);
					if (num > maxnum){
						if (num > mmx){
							mmx = num;
							cv::ellipse(image2, cv::Point(el.get_x(), el.get_y()), cv::Size(el.get_a(), el.get_b()), el.get_th(), 0, 360, cv::Scalar(255, 30, 0), 1, 1);
							cv::circle(image2, cv::Point(el.get_x(), el.get_y()), 1, cv::Scalar(50, 0, 50), -1, CV_AA);
							cv::resize(image2, image2, cv::Size(), 4, 4);
							cv::imshow("image", image2);
							/*					std::string s;
							std::string p = "_";
							s = "ellipse" + el.*/


							cv::imwrite("elli.png", image2);
							cv::waitKey(60);
						}
						std::cout << "mmx" << mmx << std::endl;
						src->copyTo(image2);
						std::cout << num << std::endl;
						el.set_count(num);
						v.push_back(el);
						//cv::ellipse(image2, cv::Point(el.get_x(), el.get_y()), Size(el.get_a(), el.get_b()), el.get_th(), 0, 360, Scalar(255, 30, 0), 1, 1);
						//cv::circle(image2, cv::Point(el.get_x(), el.get_y()), 1, Scalar(50, 0, 50), -1, CV_AA);

					}
					num = 0;
					////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				}
			}
		}
	}
}

#include "ImgProcess.h"


cv::Mat element = cv::Mat::ones(3,3,CV_8UC1); //dilate�����ɕK�v�ȍs��

int vote_k[181] = { 0 };
int r_array[500] = { 0 };

void OnImageProcessing(cv::Mat3b* src, cv::Mat1b* bin_img, cv::Mat1b* bin_origin){
	//cv::Mat s = cv::imread("../snap3.png", 1);
	//cv::Mat r;
	//s.copyTo(origin);
	//cv::resize(s, s, cv::Size(), 0.05, 0.05);
	//s.copyTo(test_img);

	cv::Mat1b Bin;


	/*************�G�b�W�_�̊i�[���s��**************/
	std::vector<cv::Point2i> EdgeArray;
	storeEdges(EdgeArray, bin_img); //�G�b�W�_���i�[
	std::cout <<"Total Edge number : "<< EdgeArray.size() << std::endl; //�i�[���ꂽ�G�b�W�_����\��



	/*****����G�b�W�_�𒍖ړ_�Ƃ��Ă��̎���13x13�̂�����2�y�AEi,Ej��ݒ肵�A
										�����Œ��ڃG�b�W�_�̐ڐ��̕����������߂�****************/
	std::vector<Tangent> TV;

	std::cout << "start deriveTangent " << std::endl;

	std::vector<Tangent> tanVec = deriveTangent(TV, EdgeArray, bin_img);

	std::cout << "end deriveTangent " << std::endl;





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


void storeEdges(std::vector<cv::Point2i> &EdgeArray, cv::Mat1b * bin_img){
	
	for (int j = 0; j < bin_img->rows; j++){
		unsigned char *ptr = bin_img->ptr<unsigned char>(j);
		for (int i = 0; i < bin_img->cols; i++){
			std::cout << " " << (int)ptr[i] / 250;  //�R���\�[���\���̍ۂ�
			//�s�N�Z���l��0���傫���ꍇ && �摜��6�s�N�Z�������ɓ����Ƀs�N�Z����
			//���݂��Ă���΂��̃s�N�Z�����i�[����
			if ((int)ptr[i] > 0 &&    
				6 < j && j < bin_img->rows - 6 &&
				6 < i && i < bin_img->cols - 6){
				cv::Point2i edgeCoordinate_xy(i, j);
				EdgeArray.push_back(edgeCoordinate_xy);
			}
		}
		std::cout << std::endl;//�\���̂��߂̉��s
	}
	

}


//�ڐ��̏����擾
std::vector<Tangent> deriveTangent(std::vector<Tangent> &Tan, std::vector<cv::Point2i> &E, cv::Mat1b *img){
	std::vector<int> rad; //�~�̔��a�i�[�p�x�N�^�z��
	std::vector<cv::Point2i>::iterator it;
	std::vector<cv::Point2i> EdgeSet;
	double s = 0;
	cv::Point2i AE, AB, e1, e2;
	cv::Point2i ate;


	//�i�[���ꂽ���ׂẴG�b�W�ɑ΂��ď������s��
	for (auto it : E){
		int vote_k[361] = { 0 };
		int r_array[500] = { 0 };

		//���ڃG�b�W�_�ɑ΂���13x13�̗̈������
		for (int j = it.y - 6; j < it.y + 7; j++){
			unsigned char *ptr = img->ptr<unsigned char>(j);
			for (int i = it.x - 6; i < it.x + 7; i++){
				//@comment ���S�Ɠ����ꍇ�A�摜�̗̈悩��O���ꍇ��pass
				if ((j != it.y && i != it.x) &&
					(j > 0 && j < img->rows) &&
					(i > 0 && i < img->cols)
					){
					if ((int)ptr[i] > 0){
						cv::Point2i near_edge = cv::Point2i(it.x, it.y);
						//ce.x = it.x;
						//ce.y = it.y;
						//@comment 8�ߖT�ɑ��݂���G�b�W�͏Ȃ��B
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
		derive_Theta_R(EdgeSet, attentionEdge, vote_k, rad); //�~�̔��a�Ɗp�x�𓊕[�ɂ���ċ��߂�


		int max_count = 1;
		int index = 0;
		int start = 0;

		//�����Ƃ��J�E���g���傫���p�x�����̉~�̊p�x�Ƃ���
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
		//�ł��J�E���g���傫�����a�����̉~�̔��a�Ƃ���
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
		//Tangent�N���X�ɒ��ڃG�b�W�̏���o�^����
		if (index != 0 && r_index != 0){
			Tangent t = Tangent(attentionEdge, index - 180, r_index);
			Tan.push_back(t);
		}
		EdgeSet.clear();
	}
	return Tan;
}


//@comment 8�ߖT�̑��݂��`�F�b�N
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

//����G�b�W�_�ł̃Ƃ�R�𓊕[�ɂ���ċ��߂�
void derive_Theta_R(std::vector<cv::Point2i> &edges, cv::Point2i &ate, int *vote_k, std::vector<int> &rad){
	std::vector<cv::Point2i>::iterator vb, ve; //�z��p�C�e���[�^
	cv::Point2i e1, e2, center = ate;
	cv::Point2i AE, AB;
	double s = 0;

	//@comment ���ׂẴG�b�W�g�ݍ��킹�ɑ΂��ď������s��
	for (auto vb : edges){
		for (auto ve : edges){
			//@comment �����G�b�W�_�͏������Ȃ�
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
				//ate(attention_edge)��Pi,Pj�̂Ȃ��p���U�O�x�ȉ��̏ꍇ�͏������s��Ȃ��B
				if (s >= 60){
					calculate_center(&center, &e1, &e2, vote_k, rad);
				}
			}
		}
	}
}


//@comment �^����ꂽ�R�_����~�̒��S�𐄒肷��
// �s�N�Z���̗ʎq�덷���l�����ăs�N�Z���̂S����ݒ肷��
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

	//3�_��ʂ�~�̒��S�����߂�
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
			C = A.inv() * B;//�s�񉉎Z�ɂ����W�擾

			it_x = C.begin<float>();
			it_y = C.end<float>() - 1;
			cpara.x = *it_x;
			cpara.y = *it_y;
			center.push_back(cpara);
			//�~�̔��a�擾
			r = sqrt(((c->x - (*it_x)) *(c->x - (*it_x))) + ((c->y - (*it_y))*(c->y - (*it_y))));
			radius.push_back(r);
		}
	}


	vote_theta(center, *c, vote);//�p�x�Ɋւ��铊�[

	vote_rad(radius, rad);// ���a�ɑ΂��铊�[
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
		if ((it.x - P.x)> 0 && (it.y - P.y)< 0){ //���ی�
			theta = -atan2(it.y - P.y, it.x - P.x) * 180 / CV_PI;
		}
		else if ((it.x - P.x)< 0 && (it.y - P.y) < 0){//���ی�
			theta = -atan2(it.y - P.y, it.x - P.x) * 180 / CV_PI;
			if (theta < 0){
			}
		}
		else if ((it.x - P.x) > 0 && (it.y - P.y) > 0){//��l�ی�
			theta = -atan2(it.y - P.y, it.x - P.x) * 180 / CV_PI;
		}
		else if ((it.x - P.x) < 0 && (it.y - P.y) > 0){//��O�ی�
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
	//vote_theta��
	for (int i = min_theta; i < max_theta; i++){
		*(vote + 180 + i) += 1;
	}
	for (int i = min_theta2; i < max_theta2; i++){
		*(vote + 180 + i) += 1;
	}
}
void vote_rad(std::vector<float> radi, std::vector<int> &r){
	//std::vector<float>::iterator so;
	std::sort(radi.begin(), radi.end()); //r���~���Ƀ\�[�g
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





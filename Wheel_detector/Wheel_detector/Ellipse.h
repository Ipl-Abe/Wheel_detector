#include <opencv/cxcore.h>

class El{
private:
	int x, y, th, a, b;
	int count;
public:
	El(int x, int y, int th, int a, int b);
	int get_x();
	int get_y();
	int get_th();
	int get_a();
	int get_b();
	int get_count();
	int set_count(int num);

};

class El_para{
private:
	cv::Point2i O;
	int a_2, b_2;
	int h, b, c;
public:
	El_para(cv::Point2i O);

	/*void set_a();
	void set_b();
	void set_h();
	void set_b();
	void set_c()
	*/;
};

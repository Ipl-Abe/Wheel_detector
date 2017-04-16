// SSK.h : SSK.DLL のメイン ヘッダー ファイル
//
#ifndef SSK_H
#define SSK_H

#pragma once

#define SUCCESS					0
#define CONNECT_ERROR			1
#define INIT_ERROR				2
#define GET_COMMAND_ERROR		-1

#define COMMAND_ERROR			1

#define SET_FUNC_ERROR			0
#define SET_FUNC_SUCCESS		1
#define FRAME_GET_ERROR			0
#define FRAME_GET_SUCCESS		1

/****** Callback Function typedef ******/
#ifndef nudecl

#ifdef WIN32
  #define nudecl __stdcall
#else
  #define nudecl
#endif

typedef INT (nudecl *ImagecallbackType)(INT c, INT no, INT nSize, UCHAR* pData, USHORT Width, USHORT Height, USHORT Format, SYSTEMTIME Time);
#endif /* nuRetrn */
/****** Callback Function typedef ******/

/****** JPEG Image Struct ******/
typedef struct {
	BYTE*		StreamImageBytes;
	long		StreamLength;
    long		ImageWidth;
    long		ImageHeight;
	SYSTEMTIME	Time;
	int			no;
} ImageStruct;
/****** JPEG Image Struct ******/

/****** DISPLAY Mode ******/
#define RING			0
#define PANORAMA1		1
#define PANORAMA1_FLIP	2
#define PANORAMA2		3
#define PANORAMA2_FLIP	4
#define PANORAMA4		5
#define PANORAMA4_FLIP	6
#define ORTHOGONALIZE	11
#define LARGEMODE		9
/****** DISPLAY Mode ******/

/****** MOVE TYPE ******/
#define TYPE_AVI	0//AVI形式
#define TYPE_MJPEG	1//JPEG形式
/****** MOVE TYPE ******/

/****** DLL Function ******/
#ifdef __cplusplus
extern "C"{
#endif
int nm30_init(void);
int nm30_get_number(void);
char* nm30_get_id(int c);
int nm30_param_save(int c);
int nm30_set_panorama_mode(int c, int mode);
int nm30_get_panorama_mode(int c);
int nm30_set_pan(int c, float pan);
float nm30_get_pan(int c);
int nm30_set_tilt(int c, float tilt);
float nm30_get_tilt(int c);
int nm30_set_zoom(int c, float zoom);
float nm30_get_zoom(int c);
int nm30_set_roll(int c, float roll);
float nm30_get_roll(int c);
int nm30_enable_information_display(int c);
int nm30_disable_information_display(int c);
int nm30_enable_overlay_display(int c);
int nm30_disable_overlay_display(int c);
int nm30_set_capture_size(int c, int width, int height);
int nm30_get_capture_width(int c);
int nm30_get_capture_height(int c);
int nm30_set_capture_fps(int c, float fps);
float nm30_get_capture_fps(int c);
float nm30_get_actual_fps(int c);
int nm30_jpeg_set(int c, int jpegcomp);
int nm30_jpeg_get(int c);
int nm30_set_autopan_speed(int c, int speed);
int nm30_get_autopan_speed(int c);
int nm30_set_flip_screen(int c, int mode);
int nm30_get_flip_screen(int c);
int nm30_set_sharpness(int c, int filter);
int nm30_get_sharpness(int c);
int nm30_set_exposure_time(int c, int time);
int nm30_get_exposure_time(int c);
int nm30_set_gain(int c, int gain);
int nm30_get_gain(int c);
int nm30_enable_auto_exposure(int c);
int nm30_disable_auto_exposure(int c);
int nm30_start_capture(int c);
int nm30_grab_frame(int c, int mode);
ImageStruct* nm30_retrieve_frame(int c);
int nm30_stop_capture(int c);
int nm30_register_callback(int c, ImagecallbackType func);
int nm30_start_movie_capture(int c, char* filepath, int filetype);
int nm30_stop_movie_capture(int c);
int nm30_disconnect(void);

// CamSync
int nm30_sync_set_sync_mode(int c1, int c2);
int nm30_sync_set_normal_mode(void);
int nm30_sync_get_mode(int c);
int nm30_sync_trigger(int c, int mode);
int nm30_sync_get_error(void);

#ifdef __cplusplus
}
#endif
#endif /* SSK_H */
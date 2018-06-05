#ifndef US_IMG_H
#define US_IMG_H

/**
 * @brief struct for sonography data
 * 
 */
typedef struct
{
    int line_cnt;   /* number of lines scanned */
    int spl;        /* number of sample points per line */
    float depth;
    float angle;
    float radius;
    float *pixels;  /* [line_cnt][spl] */
} USImage;

typedef struct
{
    int width;
    int height;
    unsigned char *pixels;
}GrayImage;


USImage *usi_load(const char *filename);

int gi_write_png(const char* filename, const GrayImage *gi);
int gi_write_bmp(const char* filename, const GrayImage *gi);

void usi_free(USImage *usi);

void gi_free(GrayImage *gi);

/* no interpolation or resize */
GrayImage *usi2gi(const USImage *usi);

/* cubic interpolation */
GrayImage *usi_itp_nearest(const USImage *usi);
GrayImage *usi_itp_bilinear(const USImage *usi);
GrayImage *usi_itp_bicubic_catmull_rom_spline(const USImage *usi);
/**
 * @brief Linear interpolation in each column and cubic interpolation in each row 
 * 
 * @param usi 
 * @return GrayImage* 
 */
GrayImage *usi_itp_col_linear_row_cubic(const USImage *usi);
GrayImage *usi_itp_bicubic_conv(const USImage *usi);
GrayImage *usi_itp_bicubic(const USImage *usi);

/* */
int usi_itp_png(const USImage *usi, const char *fileout, GrayImage *(*usi_itp)(const USImage *));

#ifdef USE_OPENCL
#include "ocl_utils.h"

int usi_ocl_setup(OCLResrc *ocl_resrc);
int usi_ocl_release(OCLResrc *ocl_resrc);

GrayImage *usi_itp_nearest_ocl(const USImage *usi, const OCLResrc *ocl_resrc);
/**
 * @brief Nearest interpolation using OpenCL with image API
 * 
 * @param usi 
 * @param ocl_resrc 
 * @return GrayImage* 
 */
GrayImage *usi_itp_nearest_ocl_img(const USImage *usi, const OCLResrc *ocl_resrc);
GrayImage *usi_itp_bilinear_ocl(const USImage *usi, const OCLResrc *ocl_resrc);
GrayImage *usi_itp_bicubic_ocl(const USImage *usi, const OCLResrc *ocl_resrc);
GrayImage *usi_itp_bicubic_catmull_rom_spline_ocl(const USImage *usi, const OCLResrc *ocl_resrc);
/**
 * @brief Linear interpolation in each column and cubic interpolation in each row, using OpenCL 
 * 
 * @param usi 
 * @return GrayImage* 
 */
GrayImage *usi_itp_col_linear_row_cubic_ocl(const USImage *usi, const OCLResrc *ocl_resrc);
GrayImage *usi_itp_bicubic_conv_ocl(const USImage *usi, const OCLResrc *ocl_resrc);

int usi_itp_png_ocl(const USImage *usi, const char *fileout, GrayImage *(*usi_itp)(const USImage *, const OCLResrc *), const OCLResrc *ocl_resrc);
#endif

#endif /* US_IMG_H*/
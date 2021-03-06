#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include "config.h"
#include "us_img.h"
#include "debug.h"

#ifdef TIMING
#include <time.h>

#define CPU_TIMING_START()        \
    clock_t _cpu_start, _cpu_end; \
    _cpu_start = clock();

#define CPU_TIMING_END()                                                  \
    _cpu_end = clock();                                                   \
    double _cpu_time = 1000.0 * (_cpu_end - _cpu_start) / CLOCKS_PER_SEC; \
    printf("%s: CPU time used: %.2f ms\n", __func__, _cpu_time);

#define WALL_TIMING_START()                 \
    struct timespec _wall_start, _wall_end; \
    timespec_get(&_wall_start, TIME_UTC);

#define WALL_TIMING_END()                                                                                                    \
    timespec_get(&_wall_end, TIME_UTC);                                                                                      \
    double _wall_time = 1000.0 * (_wall_end.tv_sec - _wall_start.tv_sec) + 1e-6 * (_wall_end.tv_nsec - _wall_start.tv_nsec); \
    printf("%s: Wall time used: %.2f ms\n", __func__, _wall_time);

#endif /* TIMING */

#ifdef DEBUG
#define assert_pixel(pixel)                        \
    if (!(pixel >= 0.0 && pixel <= 255.0))         \
    {                                              \
        DBG_PRINT(#pixel " invalid: %f\n", pixel); \
        exit(EXIT_FAILURE);                        \
    }
#else
#define assert_pixel(pixel)
#endif

USImage *usi_load(const char *filename)
{
    FILE *fd;
    /* return value */
    USImage *usi = NULL;

    /* MUST use binary mode ("b") on Windows */
    if (!(fd = fopen(filename, "rb")))
    {
        /* fd == NULL, no need to close
         * errno has been set
         */
#ifdef DEBUG
        perror("fopen()");
#endif
        return NULL;
    }

    usi = malloc(sizeof(USImage));
    if (5 != fscanf(fd, "%d %d %f %f %f\n",
                    &(usi->line_cnt), &(usi->spl),
                    &(usi->depth), &(usi->angle),
                    &(usi->radius)))
    {
        DBG_PRINT("Failed to read header.\n");
        goto READ_HEADER_FAILURE;
    }

    if (usi->line_cnt < 0 || usi->depth < 0 || usi->depth < 0.0 || usi->angle < 0.0 || usi->radius < 0.0)
    {
        DBG_PRINT("Invalid header.\n");
        goto READ_HEADER_FAILURE;
    }
    // DBG_PRINT("width\theight\tdepth\tangle\tradius\n");
    // DBG_PRINT("%d\t%d\t%.3f\t%.3f\t%.3f\n",
    //             usi->line_cnt, usi->spl, usi->depth, usi->angle, usi->radius);

    size_t pixel_cnt = (size_t)usi->line_cnt * (size_t)usi->spl;
    /* pixels[spl][line_cnt] */
    usi->pixels = malloc(pixel_cnt * sizeof(float));
    /* pixels_tmp[line_cnt][spl] */
    float *pixels_tmp = malloc(pixel_cnt * sizeof(float));
    size_t readpixel_cnt = fread(pixels_tmp, sizeof(float), pixel_cnt, fd);
    // DBG_PRINT("pixel count: %lu\n", pixel_cnt);
    // DBG_PRINT("read pixel count: %lu\n", readpixel_cnt);
    if (pixel_cnt > readpixel_cnt)
    {
        DBG_PRINT("Failed to read pixels.\n");
        goto READ_PIXEL_FAILURE;
    }

    for (int i = 0; i < usi->line_cnt; ++i)
    {
        for (int j = 0; j < usi->spl; ++j)
        {
            usi->pixels[j * usi->line_cnt + i] = pixels_tmp[i * usi->spl + j];
        }
    }
    /* all pixels read */
    free(pixels_tmp);
    goto SUCCESS_RETURN;

READ_PIXEL_FAILURE:
    free(pixels_tmp);
    free(usi->pixels);
    usi->pixels = NULL;
READ_HEADER_FAILURE:
    free(usi);
    usi = NULL; /* sets return value to NULL */
    if (ferror(fd))
    {
#ifdef DEBUG
        perror("File I/O error");
#endif
    }
SUCCESS_RETURN:
    fclose(fd);
    return usi;
}

void usi_free(USImage *usi)
{
    free(usi->pixels);
    free(usi);
}

void gi_free(GrayImage *gi)
{
    free(gi->pixels);
    free(gi);
}

GrayImage *usi2gi(const USImage *usi)
{
    GrayImage *gi = malloc(sizeof(GrayImage));
    gi->height = usi->spl;
    gi->width = usi->line_cnt;
    gi->pixels = malloc(gi->width * gi->height * sizeof(char));
    for (int i = 0; i < gi->height; ++i)
    {
        for (int j = 0; j < gi->width; ++j)
        {
            gi->pixels[i * gi->width + j] = (unsigned char)roundf(usi->pixels[i * usi->line_cnt + j]);
        }
    }
    return gi;
}

typedef enum
{
    BORDER_REFLECT
} border_t;
/**
 * @brief 
 * 
 * @param usi 
 * @param left column index of the top-left pixel of the patch
 * @param top row index of the top-left pixel of the patch
 * @param patch 
 * @param w 
 * @param h 
 * @param border_type 
 * @return int 
 * 
 * @attention currently, border width mustn't exceed 1, and patch mustn't be larger than usi
 */
static inline int usi_get_patch(const USImage *usi, int left, int top, float *patch, int w, int h, border_t border_type)
{
    // size limit
    assert(w > 0 && h > 0 && w <= usi->line_cnt && h <= usi->spl);
    // border limit
    assert(left >= -1 && top >= -1 && left + w - 1 <= usi->line_cnt && top + h - 1 <= usi->spl);
    switch (border_type)
    {
        int x, y;
    case BORDER_REFLECT:
        for (int i = 0; i < h; ++i)
        {
            y = top + i;
            y = y < 0 ? -y - 1 : y;
            y = y > usi->spl - 1 ? 2 * usi->spl - y - 1 : y;
            for (int j = 0; j < w; ++j)
            {
                x = left + j;
                x = x < 0 ? -x - 1 : x;
                x = x > usi->line_cnt - 1 ? 2 * usi->line_cnt - x - 1 : x;
                patch[i * w + j] = usi->pixels[y * usi->line_cnt + x];
            }
        }
        break;
    default:
        return -1;
    }
    return 0;
}

static inline float itp_linear(float y1, float y2, float dx)
{
    return (1.0f - dx) * y1 + dx * y2;
}

/**
 * @brief 
 * 
 * @param y f(0) = y[1], f(1) = y[2], f'(0) = (y[2]-y[0])/2, f'(1) = (y[3]-y[1])/2
 * @param dx 
 * @return float 
 */
static inline float itp_catmull_rom_spline(float y[4], float dx)
{
    float a[4];
    a[3] = -0.5f * y[0] + 1.5f * y[1] - 1.5f * y[2] + 0.5f * y[3];
    a[2] = y[0] - 2.5f * y[1] + 2.0f * y[2] - 0.5f * y[3];
    a[1] = -0.5f * y[0] + 0.5f * y[2];
    a[0] = y[1];
    return fmaf(fmaf(fmaf(a[3], dx, a[2]), dx, a[1]), dx, a[0]);
}

GrayImage *usi_itp_nearest(const USImage *usi)
{
#ifdef TIMING
    CPU_TIMING_START();
#endif
    float s_interval = usi->depth / (usi->spl - 1);          /* sampling interval */
    float a_interval = 2 * usi->angle / (usi->line_cnt - 1); /* angle interval */
    float R = usi->radius + usi->depth;
    float real_h = R - usi->radius * cosf(usi->angle);
    float real_w = 2 * R * sinf(usi->angle);
    float center2top = usi->radius * cosf(usi->angle);

    GrayImage *gi = malloc(sizeof(GrayImage));
    gi->height = (int)ceilf(usi->spl * real_h / usi->depth); /* uses the same sampling interval */
    gi->width = (int)ceilf(gi->height * real_w / real_h);
    /* allocates memory and sets all bits to 0 */
    gi->pixels = calloc(gi->width * gi->height, sizeof(char));

    /* forward map */
    // for(int i = 0; i < usi->spl; ++i)
    // {
    //     float real_rho = usi->radius + i * usi->depth/(usi->spl-1);
    //     for(int j = 0; j < usi->line_cnt; ++j)
    //     {
    //         float theta = -usi->angle + j * (2*usi->angle/(usi->line_cnt-1));
    //         /* suppose that the first sample point is on the real edge */
    //         float real_x = real_w/2 + real_rho * sinf(theta);
    //         float real_y = real_rho * cosf(theta) - center2top;
    //     }
    // }

    /* backward map */
    for (int i = 0; i < gi->height; ++i)
    {
        /* from top to bottom */
        float real_y = center2top + i * s_interval; /* real-world distance in y direction */
        for (int j = 0; j < gi->width; ++j)
        {
            float real_x = j * s_interval - real_w / 2; /* real-world distance in y direction */
            float theta = atanf(real_x / real_y);
            float rho = sqrtf(real_x * real_x + real_y * real_y);
            if (theta >= -usi->angle && theta <= usi->angle && rho <= R && rho >= usi->radius)
            {
                int j_0 = (int)((theta + usi->angle) / a_interval);
                int i_0 = (int)((rho - usi->radius) / s_interval);
                gi->pixels[i * gi->width + j] = (unsigned char)roundf(usi->pixels[i_0 * usi->line_cnt + j_0]);
                // gi->pixels[i*gi->width + j] = (unsigned char)roundf(usi->pixels[in_j*usi->spl + in_i]);
            }
        }
    }
#ifdef TIMING
    CPU_TIMING_END();
#endif
    return gi;
}

GrayImage *usi_itp_bilinear(const USImage *usi)
{
#ifdef TIMING
    CPU_TIMING_START();
#endif
    float s_interval = usi->depth / (usi->spl - 1);          /* sampling interval */
    float a_interval = 2 * usi->angle / (usi->line_cnt - 1); /* angle interval */
    float R = usi->radius + usi->depth;
    float real_h = R - usi->radius * cosf(usi->angle);
    float real_w = 2 * R * sinf(usi->angle);
    float center2top = usi->radius * cosf(usi->angle);

    GrayImage *gi = malloc(sizeof(GrayImage));
    gi->height = (int)ceilf(usi->spl * real_h / usi->depth); /* uses the same sampling interval */
    gi->width = (int)ceilf(gi->height * real_w / real_h);
    /* allocates memory and sets all bits to 0 */
    gi->pixels = calloc(gi->width * gi->height, sizeof(char));

    /* backward map */
    for (int i = 0; i < gi->height; ++i)
    {
        /* from top to bottom */
        float real_y = center2top + i * s_interval; /* real-world distance in y direction */
        for (int j = 0; j < gi->width; ++j)
        {
            float real_x = j * s_interval - real_w / 2; /* real-world distance in x direction */
            float theta = atanf(real_x / real_y);
            float rho = sqrtf(real_x * real_x + real_y * real_y);
            if (theta >= -usi->angle && theta <= usi->angle && rho <= R && rho >= usi->radius)
            {
                float u = (theta + usi->angle) / a_interval; /* column index */
                float v = (rho - usi->radius) / s_interval;  /* row index */
                int left = (int)u;
                int right = left + 1;
                int above = (int)v;
                int below = above + 1;

                float pixel_above_left = usi->pixels[above * usi->line_cnt + left];
                float pixel_above_right = usi->pixels[above * usi->line_cnt + right];
                float pixel_below_left = usi->pixels[below * usi->line_cnt + left];
                float pixel_below_right = usi->pixels[below * usi->line_cnt + right];

                float pixel_above = (u - left) * pixel_above_right - (u - right) * pixel_above_left;
                float pixel_below = (u - left) * pixel_below_right - (u - right) * pixel_below_left;
                float pixel_bilinear = (v - above) * pixel_below - (v - below) * pixel_above;
                unsigned char pixel = (unsigned char)roundf(pixel_bilinear);

                gi->pixels[i * gi->width + j] = pixel;
            }
        }
    }
#ifdef TIMING
    CPU_TIMING_END();
#endif
    return gi;
}

/**
 * @brief Linear interpolation in each column and cubic interpolation in each row, using OpenCL 
 * 
 * @param usi 
 * @return GrayImage* 
 */
GrayImage *usi_itp_col_linear_row_cubic(const USImage *usi)
{
#ifdef TIMING
    CPU_TIMING_START();
#endif
    float s_interval = usi->depth / (usi->spl - 1);          /* sampling interval */
    float a_interval = 2 * usi->angle / (usi->line_cnt - 1); /* angle interval */
    float R = usi->radius + usi->depth;
    float real_h = R - usi->radius * cosf(usi->angle);
    float real_w = 2 * R * sinf(usi->angle);
    float center2top = usi->radius * cosf(usi->angle);

    GrayImage *gi = malloc(sizeof(GrayImage));
    gi->height = (int)ceilf(usi->spl * real_h / usi->depth); /* uses the same sampling interval */
    gi->width = (int)ceilf(gi->height * real_w / real_h);
    /* allocates memory and sets all bits to 0 */
    gi->pixels = calloc(gi->width * gi->height, sizeof(char));

    /* backward map */
    for (int i = 0; i < gi->height; ++i)
    {
        /* from top to bottom */
        float real_y = center2top + i * s_interval; /* real-world distance in y direction */
        for (int j = 0; j < gi->width; ++j)
        {
            float real_x = j * s_interval - real_w / 2; /* real-world distance in x direction */
            float theta = atanf(real_x / real_y);
            float rho = sqrtf(real_x * real_x + real_y * real_y);
            if (theta >= -usi->angle && theta <= usi->angle && rho <= R && rho >= usi->radius)
            {
                float u = (theta + usi->angle) / a_interval; /* column index */
                float v = (rho - usi->radius) / s_interval;  /* row index */
                int p = (int)u;
                int q = (int)v;
                float dx = u - p;
                float dy = v - q;

                float f_2x4[2][4];
                if (0 > usi_get_patch(usi, p - 1, q - 1, f_2x4, 4, 2, BORDER_REFLECT))
                {
                    DBG_PRINT("Failed to get a 4*2 patch at (%d, %d).\n", p - 1, q - 1);
                    gi_free(gi);
                    gi = NULL;
                    return gi;
                }

                // linear interpolation in y direction
                float lin_y[4];
                for (int n = 0; n < 4; ++n)
                    lin_y[n] = itp_linear(f_2x4[0][n], f_2x4[1][n], dy);

                // cubic interpolation in x direction
                float cub_x = itp_catmull_rom_spline(lin_y, dx);
                cub_x = cub_x < 0.0f ? 0.0f : cub_x;
                cub_x = cub_x > 255.0f ? 255.0f : cub_x;

                gi->pixels[i * gi->width + j] = (unsigned char)cub_x;
            }
        }
    }
#ifdef TIMING
    CPU_TIMING_END();
#endif
    return gi;
}

/**
 * @brief Cubic interpolation in each row first, then linear interpolation in each column
 * 
 * @param usi 
 * @return GrayImage* 
 */
GrayImage *usi_itp_row_cubic_col_linear(const USImage *usi)
{
#ifdef TIMING
    CPU_TIMING_START();
#endif
    float s_interval = usi->depth / (usi->spl - 1);          /* sampling interval */
    float a_interval = 2 * usi->angle / (usi->line_cnt - 1); /* angle interval */
    float R = usi->radius + usi->depth;
    float real_h = R - usi->radius * cosf(usi->angle);
    float real_w = 2 * R * sinf(usi->angle);
    float center2top = usi->radius * cosf(usi->angle);

    GrayImage *gi = malloc(sizeof(GrayImage));
    gi->height = (int)ceilf(usi->spl * real_h / usi->depth); /* uses the same sampling interval */
    gi->width = (int)ceilf(gi->height * real_w / real_h);
    /* allocates memory and sets all bits to 0 */
    gi->pixels = calloc(gi->width * gi->height, sizeof(char));

    /* backward map */
    for (int i = 0; i < gi->height; ++i)
    {
        /* from top to bottom */
        float real_y = center2top + i * s_interval; /* real-world distance in y direction */
        for (int j = 0; j < gi->width; ++j)
        {
            float real_x = j * s_interval - real_w / 2; /* real-world distance in x direction */
            float theta = atanf(real_x / real_y);
            float rho = sqrtf(real_x * real_x + real_y * real_y);
            if (theta >= -usi->angle && theta <= usi->angle && rho <= R && rho >= usi->radius)
            {
                float u = (theta + usi->angle) / a_interval; /* column index */
                float v = (rho - usi->radius) / s_interval;  /* row index */
                int p = (int)u;
                int q = (int)v;
                float dx = u - p;
                float dy = v - q;

                float f_2x4[2][4];
                if (0 > usi_get_patch(usi, p - 1, q - 1, f_2x4, 4, 2, BORDER_REFLECT))
                {
                    DBG_PRINT("Failed to get a 4*2 patch at (%d, %d).\n", p - 1, q - 1);
                    gi_free(gi);
                    gi = NULL;
                    return gi;
                }

                // cubic interpolation in x direction
                float cub_x[2];
                cub_x[0] = itp_catmull_rom_spline(f_2x4[0], dx);
                cub_x[1] = itp_catmull_rom_spline(f_2x4[1], dx);

                // linear interpolation in y direction
                float lin_y = itp_linear(cub_x[0], cub_x[1], dy);

                lin_y = lin_y < 0.0f ? 0.0f : lin_y;
                lin_y = lin_y > 255.0f ? 255.0f : lin_y;

                gi->pixels[i * gi->width + j] = (unsigned char)lin_y;
            }
        }
    }
#ifdef TIMING
    CPU_TIMING_END();
#endif
    return gi;
}

/**
 * @brief Cubic interpolation in each column first, then linear interpolation in each row 
 * 
 * @param usi 
 * @return GrayImage* 
 */
GrayImage *usi_itp_col_cubic_row_linear(const USImage *usi)
{
#ifdef TIMING
    CPU_TIMING_START();
#endif
    float s_interval = usi->depth / (usi->spl - 1);          /* sampling interval */
    float a_interval = 2 * usi->angle / (usi->line_cnt - 1); /* angle interval */
    float R = usi->radius + usi->depth;
    float real_h = R - usi->radius * cosf(usi->angle);
    float real_w = 2 * R * sinf(usi->angle);
    float center2top = usi->radius * cosf(usi->angle);

    GrayImage *gi = malloc(sizeof(GrayImage));
    gi->height = (int)ceilf(usi->spl * real_h / usi->depth); /* uses the same sampling interval */
    gi->width = (int)ceilf(gi->height * real_w / real_h);
    /* allocates memory and sets all bits to 0 */
    gi->pixels = calloc(gi->width * gi->height, sizeof(char));

    /* backward map */
    for (int i = 0; i < gi->height; ++i)
    {
        /* from top to bottom */
        float real_y = center2top + i * s_interval; /* real-world distance in y direction */
        for (int j = 0; j < gi->width; ++j)
        {
            float real_x = j * s_interval - real_w / 2; /* real-world distance in x direction */
            float theta = atanf(real_x / real_y);
            float rho = sqrtf(real_x * real_x + real_y * real_y);
            if (theta >= -usi->angle && theta <= usi->angle && rho <= R && rho >= usi->radius)
            {
                float u = (theta + usi->angle) / a_interval; /* column index */
                float v = (rho - usi->radius) / s_interval;  /* row index */
                int p = (int)u;
                int q = (int)v;
                float dx = u - p;
                float dy = v - q;

                float f_4x2[4][2];
                if (0 > usi_get_patch(usi, p - 1, q - 1, f_4x2, 2, 4, BORDER_REFLECT))
                {
                    DBG_PRINT("Failed to get a 2*4 patch at (%d, %d).\n", p - 1, q - 1);
                    gi_free(gi);
                    gi = NULL;
                    return gi;
                }
                // transpose f_4x2
                float traned[2][4];
                for (int m = 0; m < 4; ++m)
                    for (int n = 0; n < 2; ++n)
                        traned[n][m] = f_4x2[m][n];

                // cubic interpolation in y direction
                float cub_y[2];
                cub_y[0] = itp_catmull_rom_spline(traned[0], dy);
                cub_y[1] = itp_catmull_rom_spline(traned[1], dy);

                // linear interpolation in x direction
                float lin_x = itp_linear(cub_y[0], cub_y[1], dx);
                lin_x = lin_x < 0.0f ? 0.0f : lin_x;
                lin_x = lin_x > 255.0f ? 255.0f : lin_x;

                gi->pixels[i * gi->width + j] = (unsigned char)lin_x;
            }
        }
    }
#ifdef TIMING
    CPU_TIMING_END();
#endif
    return gi;
}

/**
 * @brief Linear interpolation in each row first, then cubic interpolation in each column
 * 
 * @param usi 
 * @return GrayImage* 
 */
GrayImage *usi_itp_row_linear_col_cubic(const USImage *usi)
{
#ifdef TIMING
    CPU_TIMING_START();
#endif
    float s_interval = usi->depth / (usi->spl - 1);          /* sampling interval */
    float a_interval = 2 * usi->angle / (usi->line_cnt - 1); /* angle interval */
    float R = usi->radius + usi->depth;
    float real_h = R - usi->radius * cosf(usi->angle);
    float real_w = 2 * R * sinf(usi->angle);
    float center2top = usi->radius * cosf(usi->angle);

    GrayImage *gi = malloc(sizeof(GrayImage));
    gi->height = (int)ceilf(usi->spl * real_h / usi->depth); /* uses the same sampling interval */
    gi->width = (int)ceilf(gi->height * real_w / real_h);
    /* allocates memory and sets all bits to 0 */
    gi->pixels = calloc(gi->width * gi->height, sizeof(char));

    /* backward map */
    for (int i = 0; i < gi->height; ++i)
    {
        /* from top to bottom */
        float real_y = center2top + i * s_interval; /* real-world distance in y direction */
        for (int j = 0; j < gi->width; ++j)
        {
            float real_x = j * s_interval - real_w / 2; /* real-world distance in x direction */
            float theta = atanf(real_x / real_y);
            float rho = sqrtf(real_x * real_x + real_y * real_y);
            if (theta >= -usi->angle && theta <= usi->angle && rho <= R && rho >= usi->radius)
            {
                float u = (theta + usi->angle) / a_interval; /* column index */
                float v = (rho - usi->radius) / s_interval;  /* row index */
                int p = (int)u;
                int q = (int)v;
                float dx = u - p;
                float dy = v - q;

                float f_4x2[4][2];
                if (0 > usi_get_patch(usi, p - 1, q - 1, f_4x2, 2, 4, BORDER_REFLECT))
                {
                    DBG_PRINT("Failed to get a 2*4 patch at (%d, %d).\n", p - 1, q - 1);
                    gi_free(gi);
                    gi = NULL;
                    return gi;
                }

                // linear interpolation in x direction
                float lin_x[4];
                for (int n = 0; n < 4; ++n)
                    lin_x[n] = itp_linear(f_4x2[n][0], f_4x2[n][1], dx);

                // cubic interpolation in y direction
                float cub_y = itp_catmull_rom_spline(lin_x, dy);
                cub_y = cub_y < 0.0f ? 0.0f : cub_y;
                cub_y = cub_y > 255.0f ? 255.0f : cub_y;

                gi->pixels[i * gi->width + j] = (unsigned char)cub_y;
            }
        }
    }
#ifdef TIMING
    CPU_TIMING_END();
#endif
    return gi;
}

GrayImage *usi_itp_bicubic_catmull_rom_spline(const USImage *usi)
{
#ifdef TIMING
    CPU_TIMING_START();
#endif
    float s_interval = usi->depth / (usi->spl - 1);          /* sampling interval */
    float a_interval = 2 * usi->angle / (usi->line_cnt - 1); /* angle interval */
    float R = usi->radius + usi->depth;
    float real_h = R - usi->radius * cosf(usi->angle);
    float real_w = 2 * R * sinf(usi->angle);
    float center2top = usi->radius * cosf(usi->angle);

    GrayImage *gi = malloc(sizeof(GrayImage));
    gi->height = (int)ceilf(usi->spl * real_h / usi->depth); /* uses the same sampling interval */
    gi->width = (int)ceilf(gi->height * real_w / real_h);
    /* allocates memory and sets all bits to 0 */
    gi->pixels = calloc(gi->width * gi->height, sizeof(char));

    /* backward map */
    for (int i = 0; i < gi->height; ++i)
    {
        /* from top to bottom */
        float real_y = center2top + i * s_interval; /* real-world distance in y direction */
        for (int j = 0; j < gi->width; ++j)
        {
            float real_x = j * s_interval - real_w / 2; /* real-world distance in x direction */
            float theta = atanf(real_x / real_y);
            float rho = sqrtf(real_x * real_x + real_y * real_y);
            if (theta >= -usi->angle && theta <= usi->angle && rho <= R && rho >= usi->radius)
            {
                float u = (theta + usi->angle) / a_interval; /* column index */
                float v = (rho - usi->radius) / s_interval;  /* row index */
                int p = (int)u;
                int q = (int)v;
                float dx = u - p;
                float dy = v - q;

                float f_4x4[4][4];
                if (0 > usi_get_patch(usi, p - 1, q - 1, f_4x4, 4, 4, BORDER_REFLECT))
                {
                    DBG_PRINT("Failed to get a 4*4 patch at (%d, %d).\n", p - 1, q - 1);
                    gi_free(gi);
                    gi = NULL;
                    return gi;
                }

                // Catmull-Rom spline interpolation in x direction
                float cub_x[4];
                for (int n = 0; n < 4; ++n)
                    cub_x[n] = itp_catmull_rom_spline(f_4x4[n], dx);

                // Catmull-Rom spline interpolation in y direction
                float cub_y = itp_catmull_rom_spline(cub_x, dy);
                cub_y = cub_y < 0.0f ? 0.0f : cub_y;
                cub_y = cub_y > 255.0f ? 255.0f : cub_y;

                gi->pixels[i * gi->width + j] = (unsigned char)cub_y;
            }
        }
    }
#ifdef TIMING
    CPU_TIMING_END();
#endif
    return gi;
}

static inline float dot(float a[], float b[], size_t n)
{
    float sum = 0;
    for (size_t i = 0; i < n; ++i)
        sum = fmaf(a[i], b[i], sum);
    return sum;
}
GrayImage *usi_itp_bicubic(const USImage *usi)
{
#ifdef TIMING
    CPU_TIMING_START();
#endif
    float s_interval = usi->depth / (usi->spl - 1);          /* sampling interval */
    float a_interval = 2 * usi->angle / (usi->line_cnt - 1); /* angle interval */
    float R = usi->radius + usi->depth;
    float real_h = R - usi->radius * cosf(usi->angle);
    float real_w = 2 * R * sinf(usi->angle);
    float center2top = usi->radius * cosf(usi->angle);

    GrayImage *gi = malloc(sizeof(GrayImage));
    gi->height = (int)ceilf(usi->spl * real_h / usi->depth); /* uses the same sampling interval */
    gi->width = (int)ceilf(gi->height * real_w / real_h);
    /* allocates memory and sets all bits to 0 */
    gi->pixels = calloc(gi->width * gi->height, sizeof(char));

    static float M_inv[16][16] =
        {
            1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            -3, 3, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            2, -2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, -2, -1, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 1, 1, 0, 0,
            -3, 0, 3, 0, 0, 0, 0, 0, -2, 0, -1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, -3, 0, 3, 0, 0, 0, 0, 0, -2, 0, -1, 0,
            9, -9, -9, 9, 6, 3, -6, -3, 6, -6, 3, -3, 4, 2, 2, 1,
            -6, 6, 6, -6, -3, -3, 3, 3, -4, 4, -2, 2, -2, -2, -1, -1,
            2, 0, -2, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0, 0, 1, 0, 1, 0,
            -6, 6, 6, -6, -4, -2, 4, 2, -3, 3, -3, 3, -2, -1, -2, -1,
            4, -4, -4, 4, 2, 2, -2, -2, 2, -2, 2, -2, 1, 1, 1, 1};

    // #ifdef DEBUG
    //     int invalid_px_count = 0;
    //     FILE *fd = fopen("bicubic_invalid_pixels.txt", "w");
    // #endif
    /* backward map */
    for (int i = 0; i < gi->height; ++i)
    {
        /* from top to bottom */
        float real_y = center2top + i * s_interval; /* real-world distance in y direction */
        for (int j = 0; j < gi->width; ++j)
        {
            float real_x = j * s_interval - real_w / 2; /* real-world distance in x direction */
            float theta = atanf(real_x / real_y);
            float rho = sqrtf(real_x * real_x + real_y * real_y);
            if (theta >= -usi->angle && theta <= usi->angle && rho <= R && rho >= usi->radius)
            {
                float u = (theta + usi->angle) / a_interval; /* column index */
                float v = (rho - usi->radius) / s_interval;  /* row index */
                int p = (int)u;                              /* column index */
                int q = (int)v;                              /* row index */
                float dx = u - p;
                float dy = v - q;

                /* derivs = [f(0,0) f(1,0) f(0,1) f(1,1) fx(0,0) ... fy(0,0) ...fxy(0,0)...]^T */
                float derivs[16] = {0.0};
                /* coefficients of the bicubic interpolation polynomial
                 * alpha = M_inv * derivs
                 * alpha = [a00 a10 a20 a30 a01 a11 a21 a31 ... a33]^T;
                 */
                float alpha[16];

                float(*f)[2] = derivs;
                float(*fx)[2] = derivs + 4;
                float(*fy)[2] = derivs + 8;
                float(*fxy)[2] = derivs + 12;

                // f(1,0) = f[0][1], f(0,1) = f[1][0]
                f[0][0] = usi->pixels[q * usi->line_cnt + p];
                f[0][1] = usi->pixels[q * usi->line_cnt + (p + 1)];
                f[1][0] = usi->pixels[(q + 1) * usi->line_cnt + p];
                f[1][1] = usi->pixels[(q + 1) * usi->line_cnt + (p + 1)];

                /*
                 * currently, for boader pixels, all derivatives are set to 0 
                 */
                // if(q > 0 && q < usi->spl - 1 && p > 0 && p < usi->line_cnt -1)
                {
                    float f_4x4[4][4];
                    // this should be replaced by a getter with border interplation, like `copyMakeBorder` in OpenCV
                    // for(int t = 0; t < 4; ++t)
                    // memcpy(f_4x4 + t, &usi->pixels[(q+t-1)*usi->line_cnt + (p-1)], 4*sizeof (float));
                    if (0 > usi_get_patch(usi, p - 1, q - 1, f_4x4, 4, 4, BORDER_REFLECT))
                    {
                        DBG_PRINT("Failed to get a 4*4 patch at (%d, %d).\n", p - 1, q - 1);
                        gi_free(gi);
                        gi = NULL;
                        return gi;
                    }
                    for (int t = 0; t < 2; ++t)
                        for (int s = 0; s < 2; ++s)
                        {
                            // fx(s,t) = (f(s+1,t) - f(s-1,t)) / 2
                            fx[t][s] = 0.5f * (f_4x4[1 + t][1 + s + 1] - f_4x4[1 + t][1 + s - 1]);
                            // fy(s,t) = (f(s,t+1) - f(s,t-1)) / 2
                            fy[t][s] = 0.5f * (f_4x4[1 + t + 1][1 + s] - f_4x4[1 + t - 1][1 + s]);
                        }
                    for (int t = 0; t < 2; ++t)
                        for (int s = 0; s < 2; ++s)
                            // fxy(s,t) = (fx(s,t+1) - fx(s,t-1)) / 2
                            //          = (f(s+1,t+1) + f(s-1,t-1) - f(s+1,t-1) - f(s-1, t+1)) / 4
                            fy[t][s] = 0.25f * (f_4x4[1 + t + 1][1 + s + 1] + f_4x4[1 + t - 1][1 + s - 1] - f_4x4[1 + t - 1][1 + s + 1] - f_4x4[1 + t + 1][1 + s - 1]);
                }

                for (size_t n = 0; n < 16; ++n)
                    alpha[n] = dot(M_inv[n], derivs, 16);

                // aij = a[j][i]
                float(*a)[4] = alpha;
                float pixel_bicubic = 0.0f;

                float result = 0.0f;
                for (int n = 3; n >= 0; --n)
                {
                    float x_part = 0.0;
                    for (int m = 3; m >= 0; --m)
                        x_part = fmaf(x_part, dx, a[n][m]);

                    result = fmaf(result, dy, x_part);
                }

                // #ifdef DEBUG
                //     if(result > 255.0f || result < 0.0f)
                //     {
                //         invalid_px_count ++;
                //         fprintf(fd, "%.2f\n", result);
                //     }
                //     // result = fabs(result);
                //     // assert(result < 270.0);
                //     // assert_pixel(result);
                // #endif

                result = result < 0.0f ? 0.0f : result;
                result = result > 255.0f ? 255.0f : result;
                gi->pixels[i * gi->width + j] = (unsigned char)roundf(result);
            }
        }
    }

    //     DBG_PRINT("invalid pixel count: %d\n", invalid_px_count);
    // #ifdef DEBUG
    //     fclose(fd);
    // #endif

#ifdef TIMING
    CPU_TIMING_END();
#endif
    return gi;
}

/* -1 on failure 
 * 0 on success
 */
int gi_write_png(const char *filename, const GrayImage *gi)
{
    return stbi_write_png(filename, gi->width, gi->height, 1, gi->pixels, 0) == 1 ? 0 : -1;
}

/* -1 on failure 
 * 0 on success
 */
int gi_write_bmp(const char *filename, const GrayImage *gi)
{
    return stbi_write_bmp(filename, gi->width, gi->height, 1, gi->pixels) == 1 ? 0 : -1;
}

int usi_itp_png(const USImage *usi, const char *pngname, GrayImage *(*usi_itp)(const USImage *))
{
    GrayImage *gi = usi_itp(usi);
    if (gi == NULL)
    {
        DBG_PRINT("Failed to interpolate.\n");
        return -1;
    }

    if (gi_write_png(pngname, gi) < 0)
    {
        DBG_PRINT("Failed to write to png.\n ");
        gi_free(gi);
        return -1;
    }

    gi_free(gi);
    return 0;
}

#ifdef USE_OPENCL
#include "ocl_utils.h"
#define PROGRAM_SRC "interpolate.cl"

int usi_ocl_setup(OCLResrc *ocl_resrc)
{
    cl_int status = ocl_setup(PROGRAM_SRC, ocl_resrc);
    if (status != CL_SUCCESS)
    {
        DBG_PRINT("Failed to set up OpenCL resources, error code: %d\n", status);
        return -1;
    }
    return 0;
}

int usi_ocl_release(OCLResrc *ocl_resrc)
{
    cl_int status = ocl_release(ocl_resrc);
    if (status != CL_SUCCESS)
    {
        DBG_PRINT("Failed to release OpenCL resources, error code: %d\n", status);
        return -1;
    }
    return 0;
}

int usi_itp_png_ocl(const USImage *usi, const char *pngname, GrayImage *(*usi_itp_ocl)(const USImage *, const OCLResrc *), const OCLResrc *ocl_resrc)
{
    GrayImage *gi = usi_itp_ocl(usi, ocl_resrc);
    if (gi == NULL)
    {
        DBG_PRINT("Failed to interpolate.\n");
        return -1;
    }

    if (gi_write_png(pngname, gi) < 0)
    {
        DBG_PRINT("Failed to write to png.\n ");
        gi_free(gi);
        return -1;
    }

    gi_free(gi);
    return 0;
}

GrayImage *usi_itp_nearest_ocl(const USImage *usi, const OCLResrc *ocl_resrc)
{
#ifdef TIMING
    WALL_TIMING_START();
#endif
    GrayImage *gi = NULL;

    cl_int status;
    cl_kernel kernel = clCreateKernel(ocl_resrc->program, "nearest", &status);
    if (status != CL_SUCCESS)
    {
        DBG_PRINT("clCreateKernel failed, error code: %d\n", status);
        return NULL;
    }

    float s_interval = usi->depth / (usi->spl - 1);          /* sampling interval */
    float a_interval = 2 * usi->angle / (usi->line_cnt - 1); /* angle interval */
    float R = usi->radius + usi->depth;
    float real_h = R - usi->radius * cosf(usi->angle);
    float real_w = 2 * R * sinf(usi->angle);
    float center2top = usi->radius * cosf(usi->angle);

    gi = malloc(sizeof(GrayImage));
    gi->height = (int)ceilf(usi->spl * real_h / usi->depth); /* uses the same sampling interval */
    gi->width = (int)ceilf(gi->height * real_w / real_h);
    /* allocates memory and sets all bits to 0 */
    gi->pixels = calloc(gi->width * gi->height, sizeof(char));

    cl_mem buf_usi_pixels = clCreateBuffer(ocl_resrc->context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float) * usi->spl * usi->line_cnt, usi->pixels, &status);
    cl_mem buf_gi_pixels = clCreateBuffer(ocl_resrc->context, CL_MEM_WRITE_ONLY, sizeof(char) * gi->height * gi->width, NULL, &status);

    status = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&buf_usi_pixels);
    status = clSetKernelArg(kernel, 1, sizeof(float), (void *)&usi->radius);
    status = clSetKernelArg(kernel, 2, sizeof(float), (void *)&usi->angle);
    status = clSetKernelArg(kernel, 3, sizeof(int), (void *)&usi->spl);
    status = clSetKernelArg(kernel, 4, sizeof(int), (void *)&usi->line_cnt);
    status = clSetKernelArg(kernel, 5, sizeof(float), (void *)&R);
    status = clSetKernelArg(kernel, 6, sizeof(float), (void *)&center2top);
    status = clSetKernelArg(kernel, 7, sizeof(float), (void *)&real_w);
    status = clSetKernelArg(kernel, 8, sizeof(float), (void *)&s_interval);
    status = clSetKernelArg(kernel, 9, sizeof(float), (void *)&a_interval);
    status = clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *)&buf_gi_pixels);

    // Define an index space of work-items for execution.
    // A work-group size is not required, but can be used.
    size_t global_size[2] = {gi->height, gi->width}, local_size[2];
    // global_id range should equal to the length of subsums

    // number-of-groups: global_size / local_size
    // work_group_size[0] = 256;

#ifdef TIMING
    cl_queue_properties queue_props[] = {CL_QUEUE_PROPERTIES,
                                         CL_QUEUE_PROFILING_ENABLE,
                                         0};
    cl_command_queue cmd_queue = clCreateCommandQueueWithProperties(ocl_resrc->context, ocl_resrc->devices[0], queue_props, &status);
    cl_event kernel_event;
    status = clEnqueueNDRangeKernel(cmd_queue, kernel, 2, NULL,
                                    global_size, NULL, 0, NULL, &kernel_event);
#else
    cl_command_queue cmd_queue = clCreateCommandQueueWithProperties(ocl_resrc->context, ocl_resrc->devices[0], NULL, &status);
    status = clEnqueueNDRangeKernel(cmd_queue, kernel, 2, NULL,
                                    global_size, NULL, 0, NULL, NULL);
#endif /* #ifdef TIMING */

    if (status != CL_SUCCESS)
    {
        DBG_PRINT("Failed to run the kernel, error code: %d\n", status);
        gi_free(gi);
        gi = NULL;
        goto RELEASE_ALL;
    }

#ifdef TIMING
    // profiles kernel execution
    clWaitForEvents(1, &kernel_event); // wait for the kernel to be finished
    cl_ulong kernel_exec_time = ocl_get_cmd_exec_time(kernel_event);
    clReleaseEvent(kernel_event);

    printf("%s: kernel execution time: %.2f ms\n", __func__, kernel_exec_time * 1e-6);
#endif

    // Read the device output buffer to the host output array
    status = clEnqueueReadBuffer(cmd_queue, buf_gi_pixels, CL_TRUE, 0,
                                 sizeof(char) * gi->height * gi->width, gi->pixels, 0, NULL, NULL);

    if (status != CL_SUCCESS)
    {
        DBG_PRINT("Failed to run the kernel, error code: %d\n", status);
        gi_free(gi);
        gi = NULL;
        goto RELEASE_ALL;
    }

RELEASE_ALL:
    clReleaseCommandQueue(cmd_queue);
    clReleaseMemObject(buf_gi_pixels);
    clReleaseMemObject(buf_usi_pixels);
    clReleaseKernel(kernel);

#ifdef TIMING
    WALL_TIMING_END();
#endif
    return gi;
}

GrayImage *usi_itp_bilinear_ocl(const USImage *usi, const OCLResrc *ocl_resrc)
{
#ifdef TIMING
    WALL_TIMING_START();
#endif
    GrayImage *gi = NULL;

    cl_int status;

    cl_kernel kernel = clCreateKernel(ocl_resrc->program, "bilinear", &status);
    if (status != CL_SUCCESS)
    {
        DBG_PRINT("clCreateKernel failed, error code: %d\n", status);
        return NULL;
    }

    float s_interval = usi->depth / (usi->spl - 1);          /* sampling interval */
    float a_interval = 2 * usi->angle / (usi->line_cnt - 1); /* angle interval */
    float R = usi->radius + usi->depth;
    float real_h = R - usi->radius * cosf(usi->angle);
    float real_w = 2 * R * sinf(usi->angle);
    float center2top = usi->radius * cosf(usi->angle);

    gi = malloc(sizeof(GrayImage));
    gi->height = (int)ceilf(usi->spl * real_h / usi->depth); /* uses the same sampling interval */
    gi->width = (int)ceilf(gi->height * real_w / real_h);
    /* allocates memory and sets all bits to 0 */
    gi->pixels = calloc(gi->width * gi->height, sizeof(char));

    cl_mem buf_usi_pixels = clCreateBuffer(ocl_resrc->context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float) * usi->spl * usi->line_cnt, usi->pixels, &status);
    cl_mem buf_gi_pixels = clCreateBuffer(ocl_resrc->context, CL_MEM_WRITE_ONLY, sizeof(char) * gi->height * gi->width, NULL, &status);

    status = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&buf_usi_pixels);
    status = clSetKernelArg(kernel, 1, sizeof(float), (void *)&usi->radius);
    status = clSetKernelArg(kernel, 2, sizeof(float), (void *)&usi->angle);
    status = clSetKernelArg(kernel, 3, sizeof(int), (void *)&usi->spl);
    status = clSetKernelArg(kernel, 4, sizeof(int), (void *)&usi->line_cnt);
    status = clSetKernelArg(kernel, 5, sizeof(float), (void *)&R);
    status = clSetKernelArg(kernel, 6, sizeof(float), (void *)&center2top);
    status = clSetKernelArg(kernel, 7, sizeof(float), (void *)&real_w);
    status = clSetKernelArg(kernel, 8, sizeof(float), (void *)&s_interval);
    status = clSetKernelArg(kernel, 9, sizeof(float), (void *)&a_interval);
    status = clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *)&buf_gi_pixels);

    // Define an index space of work-items for execution.
    // A work-group size is not required, but can be used.
    size_t global_size[2] = {gi->height, gi->width}, local_size[2];
    // global_id range should equal to the length of subsums

    // number-of-groups: global_size / local_size
    // work_group_size[0] = 256;

#ifdef TIMING
    cl_queue_properties queue_props[] = {CL_QUEUE_PROPERTIES,
                                         CL_QUEUE_PROFILING_ENABLE,
                                         0};
    cl_command_queue cmd_queue = clCreateCommandQueueWithProperties(ocl_resrc->context, ocl_resrc->devices[0], queue_props, &status);
    cl_event kernel_event;
    status = clEnqueueNDRangeKernel(cmd_queue, kernel, 2, NULL,
                                    global_size, NULL, 0, NULL, &kernel_event);
#else
    cl_command_queue cmd_queue = clCreateCommandQueueWithProperties(ocl_resrc->context, ocl_resrc->devices[0], NULL, &status);
    status = clEnqueueNDRangeKernel(cmd_queue, kernel, 2, NULL,
                                    global_size, NULL, 0, NULL, NULL);
#endif /* #ifdef TIMING */

    if (status != CL_SUCCESS)
    {
        DBG_PRINT("Failed to run the kernel, error code: %d\n", status);
        gi_free(gi);
        gi = NULL;
        goto RELEASE_ALL;
    }

#ifdef TIMING
    // profiles kernel execution
    clWaitForEvents(1, &kernel_event); // wait for the kernel to be finished
    cl_ulong kernel_exec_time = ocl_get_cmd_exec_time(kernel_event);
    clReleaseEvent(kernel_event);

    printf("%s: kernel execution time: %.2f ms\n", __func__, kernel_exec_time * 1e-6);
#endif

    // Read the device output buffer to the host output array
    status = clEnqueueReadBuffer(cmd_queue, buf_gi_pixels, CL_TRUE, 0,
                                 sizeof(char) * gi->height * gi->width, gi->pixels, 0, NULL, NULL);

    if (status != CL_SUCCESS)
    {
        DBG_PRINT("Failed to run the kernel, error code: %d\n", status);
        gi_free(gi);
        gi = NULL;
        goto RELEASE_ALL;
    }

RELEASE_ALL:
    clReleaseCommandQueue(cmd_queue);
    clReleaseMemObject(buf_gi_pixels);
    clReleaseMemObject(buf_usi_pixels);
    clReleaseKernel(kernel);

#ifdef TIMING
    WALL_TIMING_END();
#endif

    return gi;
}

GrayImage *usi_itp_bicubic_ocl(const USImage *usi, const OCLResrc *ocl_resrc)
{
#ifdef TIMING
    WALL_TIMING_START();
#endif
    GrayImage *gi = NULL;

    cl_int status;

    cl_kernel kernel = clCreateKernel(ocl_resrc->program, "bicubic", &status);
    if (status != CL_SUCCESS)
    {
        DBG_PRINT("clCreateKernel failed, error code: %d\n", status);
        return NULL;
    }

    float s_interval = usi->depth / (usi->spl - 1);          /* sampling interval */
    float a_interval = 2 * usi->angle / (usi->line_cnt - 1); /* angle interval */
    float R = usi->radius + usi->depth;
    float real_h = R - usi->radius * cosf(usi->angle);
    float real_w = 2 * R * sinf(usi->angle);
    float center2top = usi->radius * cosf(usi->angle);

    gi = malloc(sizeof(GrayImage));
    gi->height = (int)ceilf(usi->spl * real_h / usi->depth); /* uses the same sampling interval */
    gi->width = (int)ceilf(gi->height * real_w / real_h);
    /* allocates memory and sets all bits to 0 */
    gi->pixels = calloc(gi->width * gi->height, sizeof(char));

    cl_mem buf_usi_pixels = clCreateBuffer(ocl_resrc->context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float) * usi->spl * usi->line_cnt, usi->pixels, &status);
    cl_mem buf_gi_pixels = clCreateBuffer(ocl_resrc->context, CL_MEM_WRITE_ONLY, sizeof(char) * gi->height * gi->width, NULL, &status);

    status = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&buf_usi_pixels);
    status = clSetKernelArg(kernel, 1, sizeof(float), (void *)&usi->radius);
    status = clSetKernelArg(kernel, 2, sizeof(float), (void *)&usi->angle);
    status = clSetKernelArg(kernel, 3, sizeof(int), (void *)&usi->spl);
    status = clSetKernelArg(kernel, 4, sizeof(int), (void *)&usi->line_cnt);
    status = clSetKernelArg(kernel, 5, sizeof(float), (void *)&R);
    status = clSetKernelArg(kernel, 6, sizeof(float), (void *)&center2top);
    status = clSetKernelArg(kernel, 7, sizeof(float), (void *)&real_w);
    status = clSetKernelArg(kernel, 8, sizeof(float), (void *)&s_interval);
    status = clSetKernelArg(kernel, 9, sizeof(float), (void *)&a_interval);
    status = clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *)&buf_gi_pixels);

    // Define an index space of work-items for execution.
    // A work-group size is not required, but can be used.
    size_t global_size[2] = {gi->height, gi->width}, local_size[2];
    // global_id range should equal to the length of subsums

    // number-of-groups: global_size / local_size
    // work_group_size[0] = 256;

#ifdef TIMING
    cl_queue_properties queue_props[] = {CL_QUEUE_PROPERTIES,
                                         CL_QUEUE_PROFILING_ENABLE,
                                         0};
    cl_command_queue cmd_queue = clCreateCommandQueueWithProperties(ocl_resrc->context, ocl_resrc->devices[0], queue_props, &status);
    cl_event kernel_event;
    status = clEnqueueNDRangeKernel(cmd_queue, kernel, 2, NULL,
                                    global_size, NULL, 0, NULL, &kernel_event);
#else
    cl_command_queue cmd_queue = clCreateCommandQueueWithProperties(ocl_resrc->context, ocl_resrc->devices[0], NULL, &status);
    status = clEnqueueNDRangeKernel(cmd_queue, kernel, 2, NULL,
                                    global_size, NULL, 0, NULL, NULL);
#endif /* #ifdef TIMING */

    if (status != CL_SUCCESS)
    {
        DBG_PRINT("Failed to run the kernel, error code: %d\n", status);
        gi_free(gi);
        gi = NULL;
        goto RELEASE_ALL;
    }

#ifdef TIMING
    // profiles kernel execution
    clWaitForEvents(1, &kernel_event); // wait for the kernel to be finished
    cl_ulong kernel_exec_time = ocl_get_cmd_exec_time(kernel_event);
    clReleaseEvent(kernel_event);

    printf("%s: kernel execution time: %.2f ms\n", __func__, kernel_exec_time * 1e-6);
#endif

    // Read the device output buffer to the host output array
    status = clEnqueueReadBuffer(cmd_queue, buf_gi_pixels, CL_TRUE, 0,
                                 sizeof(char) * gi->height * gi->width, gi->pixels, 0, NULL, NULL);

    if (status != CL_SUCCESS)
    {
        DBG_PRINT("Failed to run the kernel, error code: %d\n", status);
        gi_free(gi);
        gi = NULL;
        goto RELEASE_ALL;
    }

RELEASE_ALL:
    clReleaseCommandQueue(cmd_queue);
    clReleaseMemObject(buf_gi_pixels);
    clReleaseMemObject(buf_usi_pixels);
    clReleaseKernel(kernel);

#ifdef TIMING
    WALL_TIMING_END();
#endif

    return gi;
}

GrayImage *usi_itp_col_linear_row_cubic_ocl(const USImage *usi, const OCLResrc *ocl_resrc)
{
#ifdef TIMING
    WALL_TIMING_START();
#endif
    GrayImage *gi = NULL;

    cl_int status;

    cl_kernel kernel = clCreateKernel(ocl_resrc->program, "col_linear_row_cubic", &status);
    if (status != CL_SUCCESS)
    {
        DBG_PRINT("clCreateKernel failed, error code: %d\n", status);
        return NULL;
    }

    float s_interval = usi->depth / (usi->spl - 1);          /* sampling interval */
    float a_interval = 2 * usi->angle / (usi->line_cnt - 1); /* angle interval */
    float R = usi->radius + usi->depth;
    float real_h = R - usi->radius * cosf(usi->angle);
    float real_w = 2 * R * sinf(usi->angle);
    float center2top = usi->radius * cosf(usi->angle);

    gi = malloc(sizeof(GrayImage));
    gi->height = (int)ceilf(usi->spl * real_h / usi->depth); /* uses the same sampling interval */
    gi->width = (int)ceilf(gi->height * real_w / real_h);
    /* allocates memory and sets all bits to 0 */
    gi->pixels = calloc(gi->width * gi->height, sizeof(char));

    cl_mem buf_usi_pixels = clCreateBuffer(ocl_resrc->context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float) * usi->spl * usi->line_cnt, usi->pixels, &status);
    cl_mem buf_gi_pixels = clCreateBuffer(ocl_resrc->context, CL_MEM_WRITE_ONLY, sizeof(char) * gi->height * gi->width, NULL, &status);

    status = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&buf_usi_pixels);
    status = clSetKernelArg(kernel, 1, sizeof(float), (void *)&usi->radius);
    status = clSetKernelArg(kernel, 2, sizeof(float), (void *)&usi->angle);
    status = clSetKernelArg(kernel, 3, sizeof(int), (void *)&usi->spl);
    status = clSetKernelArg(kernel, 4, sizeof(int), (void *)&usi->line_cnt);
    status = clSetKernelArg(kernel, 5, sizeof(float), (void *)&R);
    status = clSetKernelArg(kernel, 6, sizeof(float), (void *)&center2top);
    status = clSetKernelArg(kernel, 7, sizeof(float), (void *)&real_w);
    status = clSetKernelArg(kernel, 8, sizeof(float), (void *)&s_interval);
    status = clSetKernelArg(kernel, 9, sizeof(float), (void *)&a_interval);
    status = clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *)&buf_gi_pixels);

    // Define an index space of work-items for execution.
    // A work-group size is not required, but can be used.
    size_t global_size[2] = {gi->height, gi->width}, local_size[2];
    // global_id range should equal to the length of subsums

    // number-of-groups: global_size / local_size
    // work_group_size[0] = 256;

#ifdef TIMING
    cl_queue_properties queue_props[] = {CL_QUEUE_PROPERTIES,
                                         CL_QUEUE_PROFILING_ENABLE,
                                         0};
    cl_command_queue cmd_queue = clCreateCommandQueueWithProperties(ocl_resrc->context, ocl_resrc->devices[0], queue_props, &status);
    cl_event kernel_event;
    status = clEnqueueNDRangeKernel(cmd_queue, kernel, 2, NULL,
                                    global_size, NULL, 0, NULL, &kernel_event);
#else
    cl_command_queue cmd_queue = clCreateCommandQueueWithProperties(ocl_resrc->context, ocl_resrc->devices[0], NULL, &status);
    status = clEnqueueNDRangeKernel(cmd_queue, kernel, 2, NULL,
                                    global_size, NULL, 0, NULL, NULL);
#endif /* #ifdef TIMING */

    if (status != CL_SUCCESS)
    {
        DBG_PRINT("Failed to run the kernel, error code: %d\n", status);
        gi_free(gi);
        gi = NULL;
        goto RELEASE_ALL;
    }

#ifdef TIMING
    // profiles kernel execution
    clWaitForEvents(1, &kernel_event); // wait for the kernel to be finished
    cl_ulong kernel_exec_time = ocl_get_cmd_exec_time(kernel_event);
    clReleaseEvent(kernel_event);

    printf("%s: kernel execution time: %.2f ms\n", __func__, kernel_exec_time * 1e-6);
#endif

    // Read the device output buffer to the host output array
    status = clEnqueueReadBuffer(cmd_queue, buf_gi_pixels, CL_TRUE, 0,
                                 sizeof(char) * gi->height * gi->width, gi->pixels, 0, NULL, NULL);

    if (status != CL_SUCCESS)
    {
        DBG_PRINT("Failed to run the kernel, error code: %d\n", status);
        gi_free(gi);
        gi = NULL;
        goto RELEASE_ALL;
    }

RELEASE_ALL:
    clReleaseCommandQueue(cmd_queue);
    clReleaseMemObject(buf_gi_pixels);
    clReleaseMemObject(buf_usi_pixels);
    clReleaseKernel(kernel);

#ifdef TIMING
    WALL_TIMING_END();
#endif

    return gi;
}

GrayImage *usi_itp_bicubic_catmull_rom_spline_ocl(const USImage *usi, const OCLResrc *ocl_resrc)
{
#ifdef TIMING
    WALL_TIMING_START();
#endif
    GrayImage *gi = NULL;

    cl_int status;

    cl_kernel kernel = clCreateKernel(ocl_resrc->program, "bicubic_catmull_rom_spline", &status);
    if (status != CL_SUCCESS)
    {
        DBG_PRINT("clCreateKernel failed, error code: %d\n", status);
        return NULL;
    }

    float s_interval = usi->depth / (usi->spl - 1);          /* sampling interval */
    float a_interval = 2 * usi->angle / (usi->line_cnt - 1); /* angle interval */
    float R = usi->radius + usi->depth;
    float real_h = R - usi->radius * cosf(usi->angle);
    float real_w = 2 * R * sinf(usi->angle);
    float center2top = usi->radius * cosf(usi->angle);

    gi = malloc(sizeof(GrayImage));
    gi->height = (int)ceilf(usi->spl * real_h / usi->depth); /* uses the same sampling interval */
    gi->width = (int)ceilf(gi->height * real_w / real_h);
    /* allocates memory and sets all bits to 0 */
    gi->pixels = calloc(gi->width * gi->height, sizeof(char));

    cl_mem buf_usi_pixels = clCreateBuffer(ocl_resrc->context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(float) * usi->spl * usi->line_cnt, usi->pixels, &status);
    cl_mem buf_gi_pixels = clCreateBuffer(ocl_resrc->context, CL_MEM_WRITE_ONLY, sizeof(char) * gi->height * gi->width, NULL, &status);

    status = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&buf_usi_pixels);
    status = clSetKernelArg(kernel, 1, sizeof(float), (void *)&usi->radius);
    status = clSetKernelArg(kernel, 2, sizeof(float), (void *)&usi->angle);
    status = clSetKernelArg(kernel, 3, sizeof(int), (void *)&usi->spl);
    status = clSetKernelArg(kernel, 4, sizeof(int), (void *)&usi->line_cnt);
    status = clSetKernelArg(kernel, 5, sizeof(float), (void *)&R);
    status = clSetKernelArg(kernel, 6, sizeof(float), (void *)&center2top);
    status = clSetKernelArg(kernel, 7, sizeof(float), (void *)&real_w);
    status = clSetKernelArg(kernel, 8, sizeof(float), (void *)&s_interval);
    status = clSetKernelArg(kernel, 9, sizeof(float), (void *)&a_interval);
    status = clSetKernelArg(kernel, 10, sizeof(cl_mem), (void *)&buf_gi_pixels);

    // Define an index space of work-items for execution.
    // A work-group size is not required, but can be used.
    size_t global_size[2] = {gi->height, gi->width}, local_size[2];
    // global_id range should equal to the length of subsums

    // number-of-groups: global_size / local_size
    // work_group_size[0] = 256;

#ifdef TIMING
    cl_queue_properties queue_props[] = {CL_QUEUE_PROPERTIES,
                                         CL_QUEUE_PROFILING_ENABLE,
                                         0};
    cl_command_queue cmd_queue = clCreateCommandQueueWithProperties(ocl_resrc->context, ocl_resrc->devices[0], queue_props, &status);
    cl_event kernel_event;
    status = clEnqueueNDRangeKernel(cmd_queue, kernel, 2, NULL,
                                    global_size, NULL, 0, NULL, &kernel_event);
#else
    cl_command_queue cmd_queue = clCreateCommandQueueWithProperties(ocl_resrc->context, ocl_resrc->devices[0], NULL, &status);
    status = clEnqueueNDRangeKernel(cmd_queue, kernel, 2, NULL,
                                    global_size, NULL, 0, NULL, NULL);
#endif /* #ifdef TIMING */

    if (status != CL_SUCCESS)
    {
        DBG_PRINT("Failed to run the kernel, error code: %d\n", status);
        gi_free(gi);
        gi = NULL;
        goto RELEASE_ALL;
    }

#ifdef TIMING
    // profiles kernel execution
    clWaitForEvents(1, &kernel_event); // wait for the kernel to be finished
    cl_ulong kernel_exec_time = ocl_get_cmd_exec_time(kernel_event);
    clReleaseEvent(kernel_event);

    printf("%s: kernel execution time: %.2f ms\n", __func__, kernel_exec_time * 1e-6);
#endif

    // Read the device output buffer to the host output array
    status = clEnqueueReadBuffer(cmd_queue, buf_gi_pixels, CL_TRUE, 0,
                                 sizeof(char) * gi->height * gi->width, gi->pixels, 0, NULL, NULL);

    if (status != CL_SUCCESS)
    {
        DBG_PRINT("Failed to run the kernel, error code: %d\n", status);
        gi_free(gi);
        gi = NULL;
        goto RELEASE_ALL;
    }

RELEASE_ALL:
    clReleaseCommandQueue(cmd_queue);
    clReleaseMemObject(buf_gi_pixels);
    clReleaseMemObject(buf_usi_pixels);
    clReleaseKernel(kernel);

#ifdef TIMING
    WALL_TIMING_END();
#endif

    return gi;
}
#endif /* USE_OPENCL */
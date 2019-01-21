#include <stdio.h>
#include <stdlib.h>

#include "config.h"

#include "debug.h"
#include "us_img.h"

void print_usage(char *command)
{
    printf("Usage: %s FILE\n", command);
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        fprintf(stderr, "%s: missing filename.\n", argv[0]);
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    char *filename = argv[1];
    USImage *usi = usi_load(filename);
    if (usi == NULL)
    {
        printf("Failed to load sonography data from \"%s\".\n", filename);
        fprintf(stderr, "Failed to load sonography data from \"%s\".\n", filename);
        return EXIT_FAILURE;
    }

    printf("Original ultrasound image:\n");
    printf("Width\tHeight\tDepth\tAngle\tRadius\t\n");
    printf("%d\t%d\t%.3f\t%.3f\t%.3f\n",
           usi->line_cnt, usi->spl,
           usi->depth, usi->angle,
           usi->radius);

    int ret_val = EXIT_SUCCESS;
    GrayImage *gi_raw = usi2gi(usi);
    if (gi_write_png("raw.png", gi_raw) < 0)
    {
        ret_val = EXIT_FAILURE;
        fprintf(stderr, "Failed to write the raw image to png.\n ");
        goto WRITE_RAW_FAILURE;
    }

    int itp_status = 0;

#define CPU_ITP(FILE_NAME, FUNC, ALGO_NAME)                                  \
    do                                                                       \
    {                                                                        \
        itp_status = usi_itp_png(usi, FILE_NAME, FUNC);                      \
        if (itp_status < 0)                                                  \
        {                                                                    \
            ret_val = EXIT_FAILURE;                                          \
            fprintf(stderr, "Failed to do " ALGO_NAME " interpolation.\n "); \
            goto ITP_FAILURE;                                                \
        }                                                                    \
    } while (0)

    CPU_ITP("nearest.png", usi_itp_nearest, "nearest neighbour");
    CPU_ITP("bilinear.png", usi_itp_bilinear, "bilinear");
    CPU_ITP("col_linear_row_cubic.png", usi_itp_col_linear_row_cubic, "col-linear-row-cubic");
    // CPU_ITP("row_cubic_col_linear.png", usi_itp_row_cubic_col_linear, "row-cubic-col-linear");
    // CPU_ITP("col_cubic_row_linear.png", usi_itp_col_cubic_row_linear, "col-cubic-row-linear");
    // CPU_ITP("row_linear_col_cubic.png", usi_itp_row_linear_col_cubic, "row-linear-col-cubic");
    CPU_ITP("bicubic_catmull_rom_spline.png", usi_itp_bicubic_catmull_rom_spline, "bicubic Catmull-Rom spline");
    CPU_ITP("bicubic.png", usi_itp_bicubic, "bicubic");

#ifdef USE_OPENCL
    OCLResrc ocl_resrc;
    if (usi_ocl_setup(&ocl_resrc) < 0)
    {
        fprintf(stderr, "Failed to set up OpenCL resources.\n");
        ret_val = EXIT_FAILURE;
        goto ITP_FAILURE;
    }

    printf("Using OpenCL GPU...\n");

#define OCL_ITP(FILE_NAME, FUNC, ALGO_NAME)                                               \
    do                                                                                    \
    {                                                                                     \
        itp_status = usi_itp_png_ocl(usi, FILE_NAME, FUNC, &ocl_resrc);                   \
        if (itp_status < 0)                                                               \
        {                                                                                 \
            ret_val = EXIT_FAILURE;                                                       \
            fprintf(stderr, "Failed to do " ALGO_NAME " interpolation using OpenCL.\n "); \
            goto OCL_ITP_FAILURE;                                                         \
        }                                                                                 \
    } while (0)

    OCL_ITP("nearest_ocl.png", usi_itp_nearest_ocl, "nearest neighbour");
    OCL_ITP("bilinear_ocl.png", usi_itp_bilinear_ocl, "bilinear");
    OCL_ITP("col_linear_row_cubic_ocl.png", usi_itp_col_linear_row_cubic_ocl, "col-linear-row-cubic");
    OCL_ITP("bicubic_catmull_rom_spline_ocl.png", usi_itp_bicubic_catmull_rom_spline_ocl, "bicubic Catmull-Rom spline");
    OCL_ITP("bicubic_ocl.png", usi_itp_bicubic_ocl, "bicubic");

OCL_ITP_FAILURE:
    usi_ocl_release(&ocl_resrc);
    goto ITP_FAILURE;
#endif /* USE_OPENCL */

WRITE_RAW_FAILURE:
    gi_free(gi_raw);
ITP_FAILURE:
    usi_free(usi);
    return ret_val;
}

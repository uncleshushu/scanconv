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
    if(argc < 2)
    {
        fprintf(stderr, "%s: missing filename.\n", argv[0]);
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    char *filename =  argv[1];
    USImage *usi = usi_load(filename);
    if(usi == NULL)
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
    // if(gi_write_bmp("raw.bmp", gi_raw) < 0)
    if(gi_write_png("raw.png", gi_raw) < 0)
    {
        ret_val = EXIT_FAILURE;
        fprintf(stderr, "Failed to write the raw image to png.\n ");
        goto WRITE_RAW_FAILURE;
    }
    gi_free(gi_raw);
    gi_raw = NULL;

    int itp_status = 0;
    itp_status = usi_itp_png(usi, "nearest.png", usi_itp_nearest);
    if(itp_status < 0)
    {
        ret_val = EXIT_FAILURE;
        fprintf(stderr, "Failed to do nearest neighbour interpolation.\n ");
        goto ITP_FAILURE;
    }
    
    itp_status = usi_itp_png(usi, "bilinear.png", usi_itp_bilinear);
    if(itp_status < 0)
    {
        ret_val = EXIT_FAILURE;
        fprintf(stderr, "Failed to do bilinear interpolation.\n ");
        goto ITP_FAILURE;
    }

    itp_status = usi_itp_png(usi, "col_linear_row_cubic.png", usi_itp_col_linear_row_cubic);
    if(itp_status < 0)
    {
        ret_val = EXIT_FAILURE;
        fprintf(stderr, "Failed to do col-linear-row-cubic interpolation.\n ");
        goto ITP_FAILURE;
    }

    itp_status = usi_itp_png(usi, "bicubic_catmull_rom_spline.png", usi_itp_bicubic_catmull_rom_spline);
    if(itp_status < 0)
    {
        ret_val = EXIT_FAILURE;
        fprintf(stderr, "Failed to do bicubic Catmull-Rom spline interpolation.\n ");
        goto ITP_FAILURE;
    }

    itp_status = usi_itp_png(usi, "bicubic.png", usi_itp_bicubic);
    if(itp_status < 0)
    {
        ret_val = EXIT_FAILURE;
        fprintf(stderr, "Failed to do bicubic interpolation.\n ");
        goto ITP_FAILURE;
    }


#ifdef USE_OPENCL
    OCLResrc ocl_resrc;
    if(usi_ocl_setup(&ocl_resrc) < 0)
    {
        fprintf(stderr, "Failed to set up OpenCL resources.\n");
        ret_val = EXIT_FAILURE;
        goto ITP_FAILURE;
    }

    printf("Using OpenCL GPU...\n");
    
    itp_status = usi_itp_png_ocl(usi, "nearest_ocl.png", usi_itp_nearest_ocl, &ocl_resrc);
    if(itp_status < 0)
    {
        ret_val = EXIT_FAILURE;
        fprintf(stderr, "Failed to do nearest neighbour interpolation using OpenCL.\n ");
        goto OCL_ITP_FAILURE;
    }

    itp_status = usi_itp_png_ocl(usi, "bilinear_ocl.png", usi_itp_bilinear_ocl, &ocl_resrc);
    if(itp_status < 0)
    {
        ret_val = EXIT_FAILURE;
        fprintf(stderr, "Failed to do bilinear interpolation using OpenCL.\n ");
        goto OCL_ITP_FAILURE;
    }
    
    itp_status = usi_itp_png_ocl(usi, "col_linear_row_cubic_ocl.png", usi_itp_col_linear_row_cubic_ocl, &ocl_resrc);
    if(itp_status < 0)
    {
        ret_val = EXIT_FAILURE;
        fprintf(stderr, "Failed to do col-linear-row-cubic interpolation using OpenCL.\n ");
        goto OCL_ITP_FAILURE;
    }

    itp_status = usi_itp_png_ocl(usi, "bicubic_catmull_rom_spline_ocl.png", usi_itp_bicubic_catmull_rom_spline_ocl, &ocl_resrc);
    if(itp_status < 0)
    {
        ret_val = EXIT_FAILURE;
        fprintf(stderr, "Failed to do bicubic Catmull-Rom spline interpolation using OpenCL.\n ");
        goto OCL_ITP_FAILURE;
    }

    itp_status = usi_itp_png_ocl(usi, "bicubic_ocl.png", usi_itp_bicubic_ocl, &ocl_resrc);
    if(itp_status < 0)
    {
        ret_val = EXIT_FAILURE;
        fprintf(stderr, "Failed to do bicubic interpolation using OpenCL.\n ");
        goto OCL_ITP_FAILURE;
    }
    

OCL_ITP_FAILURE:
    usi_ocl_release(&ocl_resrc);
    goto ITP_FAILURE;
#endif /* USE_OPENCL */

    // if(gi_itp == NULL)
    // {
    //     DBG_PRINT("Failed to interpolate.\n");
    //     ret_val = EXIT_FAILURE;
    //     goto ITP_FAILURE;
    // }

    // // if(gi_write_bmp("interpolated.bmp", gi_itp) < 0)
    // if(gi_write_png("nearest.png", gi_itp_nearest) < 0)
    // {
    //     ret_val = EXIT_FAILURE;
    //     DBG_PRINT("Failed to write the interpolated image to png.\n ");
    // }
    // else
    //     ret_val = EXIT_SUCCESS;
    
    // gi_free(gi_itp_nearest);

WRITE_RAW_FAILURE:
    gi_free(gi_raw);
ITP_FAILURE:
    usi_free(usi);
    return ret_val;
}



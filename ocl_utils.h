#ifndef OCL_UTILS_H
#define OCL_UTILS_H

#include <ctype.h>
#include <string.h>

/* OpenCL headers all-in-one */
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif

#include "debug.h"

#define INTEL_VENDOR_ID 0x8086
#define AMD_VENDOR_ID 0x1002
#define NVIDIA_VENDOR_ID 0x10DE

typedef struct
{
    cl_device_id *devices;
    cl_context context;
    cl_program program;
} OCLResrc;

static cl_program ocl_build_from_file(char *filename, cl_context context,
                                      cl_device_id *devices, cl_uint device_cnt,
                                      cl_int *status);

static cl_int ocl_get_best_gpus(cl_uint num_entries, cl_device_id *best_gpus, cl_uint *num_best_gpus);
static cl_ulong ocl_get_cmd_exec_time(cl_event event);
static cl_int ocl_setup(char *filename, OCLResrc *ocl_resrc);
static cl_int ocl_free(OCLResrc *ocl_resrc);

/**
 * @brief Get the best GPUs
 * 
 * @param num_entries number of gpus to add to the `best_gpus` list
 * @param best_gpus list of best gpus, ignored if `NULL` 
 * @param num_best_gpus number of best gpus available, ignored if `NULL`
 * @return cl_int the status code
 */
static cl_int ocl_get_best_gpus(cl_uint num_entries, cl_device_id *best_gpus, cl_uint *num_best_gpus)
{
    cl_int status;
    // Retrieve the number of platforms
    cl_uint platform_cnt = 0;
    status = clGetPlatformIDs(0, NULL, &platform_cnt);
    // exits if no platform detected
    if (status != CL_SUCCESS || platform_cnt == 0)
    {
        DBG_PRINT("No OpenCL platform detected.\n");
        return status;
    }

    cl_platform_id *platforms = NULL;
    platforms = malloc(platform_cnt * sizeof(cl_platform_id));
    // fill the list of platforms
    status = clGetPlatformIDs(platform_cnt, platforms, NULL);
    if (status != CL_SUCCESS)
    {
        DBG_PRINT("Failed to get OpenCL platform IDs.\n");
        free(platforms);
        return status;
    }

    // find the best platform
    cl_platform_id best_platform;
    for (size_t i = 0; i < platform_cnt; ++i)
    {
        best_platform = platforms[i];
        size_t platform_name_size;
        status = clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, 0, NULL, &platform_name_size);
        if (status != CL_SUCCESS)
        {
            DBG_PRINT("Failed to get the size of platform[%zu]'s name.\n", i);
            free(platforms);
            return status;
        }
        char *platform_name = malloc(platform_name_size);
        status = clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, platform_name_size, platform_name, NULL);
        if (status != CL_SUCCESS)
        {
            DBG_PRINT("Failed to get the name of platform[%zu].\n", i);
            free(platform_name);
            free(platforms);
            return status;
        }

        for (size_t i = 0; i < platform_name_size; ++i)
        {
            platform_name[i] = toupper(platform_name[i]);
        }
        if (strncmp("INTEL", platform_name, 5) != 0)
        {
            // found discrete gpu
            DBG_PRINT("Best GPU platform found: %s\n", platform_name);
            free(platform_name);
            break;
        }

        free(platform_name);

        // // Retrieve the number of gpus on platfrom[i]
        // cl_uint gpu_cnt = 0;
        // status = clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_GPU, 0, NULL, &gpu_cnt);
        // if(status != CL_SUCCESS || gpu_cnt == 0)
        // {
        //     continue;
        // }
        // // Allocate enough space for gpus
        // gpu_groups[i] = (cl_device_id*)malloc(gpu_cnt * sizeof(cl_device_id));
        // // gets device ids of platforms[i]'s gpus
        // status = clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_GPU, gpu_cnt, gpu_groups[i], NULL);
        // if(status != CL_SUCCESS)
        // {
        //     DBG_PRINT("Failed to get device ids of platform[%d]'s gpus.\n", i);
        //     continue;
        // }
        // cl_int vendor_id;
        // status = clGetDeviceInfo(gpus_groups[i][0], CL_DEVICE_VENDOR_ID, sizeof(vendor_id), &vendor_id, NULL);
    }

    if (num_best_gpus != NULL)
    {
        // Retrieve the number of gpus on best_platform
        status = clGetDeviceIDs(best_platform, CL_DEVICE_TYPE_GPU, 0, NULL, num_best_gpus);
        if (status != CL_SUCCESS || *num_best_gpus == 0)
        {
            DBG_PRINT("No OpenCL GPU detected.\n");
            free(platforms);
            return status;
        }
    }

    if (best_gpus != NULL)
    {
        status = clGetDeviceIDs(best_platform, CL_DEVICE_TYPE_GPU, num_entries, best_gpus, NULL);
        if (status != CL_SUCCESS)
        {
            DBG_PRINT("Failed to get device IDs of best gpus.\n");
            free(platforms);
            return status;
        }
    }

    free(platforms);
    return status;
}

static cl_program ocl_build_from_file(char *filename, cl_context context,
                                      cl_device_id *devices, cl_uint device_cnt,
                                      cl_int *status)
{
    /* read the source of the opencl program */
    FILE *program_fd = fopen(filename, "r");
    fseek(program_fd, 0, SEEK_END);
    size_t program_size = ftell(program_fd);
    rewind(program_fd);
    char *program_src = malloc(program_size + 1);
    fread(program_src, sizeof(char), program_size, program_fd);
    program_src[program_size] = '\0';
    fclose(program_fd);

    // creates and compiles the program from the source code strings
    // only 1 string here, null-terminated
    cl_program program = clCreateProgramWithSource(context, 1, &program_src, NULL, status);
    free(program_src);
    if (*status != CL_SUCCESS || program == NULL)
    {
        DBG_PRINT("clCreateProgramWithSource failed.\n");
        return NULL;
    }
    char *build_options = "-cl-mad-enable -cl-std=CL2.0";
    // char *build_options = "-cl-mad-enable -cl-std=CL1.2";
    *status = clBuildProgram(program, device_cnt, devices, build_options, NULL, NULL);
    if (*status != CL_SUCCESS)
    {
        DBG_PRINT("clBuildProgram failed, error code: %d\n", *status);
        cl_build_status build_status;
        char build_log[0x10000];
        for (size_t i = 0; i < device_cnt; ++i)
        {
            *status = clGetProgramBuildInfo(program, devices[i], CL_PROGRAM_BUILD_STATUS, sizeof build_status, &build_status, NULL);
            if (*status != CL_SUCCESS)
            {
                DBG_PRINT("\nFailed to get build status of device %zu\n", i);
                continue;
            }
            if (build_status == CL_BUILD_NONE)
            {
                DBG_PRINT("\nNo build operations have been performed on the specified program for device %zu\n", i);
                continue;
            }
            *status = clGetProgramBuildInfo(program, devices[i], CL_PROGRAM_BUILD_LOG, sizeof build_log, build_log, NULL);
            DBG_PRINT("\nBuild log for device %zu:\n%s\n", i, build_log);
        }
        clReleaseProgram(program);
        return NULL;
    }
    return program;
}

// return the execution time of an event in nanoseconds
static cl_ulong ocl_get_cmd_exec_time(cl_event event)
{
    cl_ulong evt_start; // event start time
    cl_ulong evt_end;   // event end time
    clGetEventProfilingInfo(event,
                            CL_PROFILING_COMMAND_START,
                            sizeof(cl_ulong),
                            &evt_start,
                            0);
    clGetEventProfilingInfo(event,
                            CL_PROFILING_COMMAND_END,
                            sizeof(cl_ulong),
                            &evt_end,
                            0);
    return evt_end - evt_start;
}

static cl_int ocl_setup(char *filename, OCLResrc *ocl_resrc)
{
    cl_int status;
    cl_uint gpu_cnt = 0;
    status = ocl_get_best_gpus(0, NULL, &gpu_cnt);
    if (status != CL_SUCCESS || gpu_cnt == 0)
    {
        DBG_PRINT("Failed to get OpenCL GPU.\n");
        return status;
    }
    ocl_resrc->devices = malloc(sizeof(cl_device_id) * gpu_cnt);
    status = ocl_get_best_gpus(gpu_cnt, ocl_resrc->devices, NULL);
    if (status != CL_SUCCESS)
    {
        DBG_PRINT("Failed to get OpenCL GPU.\n");
        goto FREE_GPUS;
    }

    ocl_resrc->context = clCreateContext(NULL, gpu_cnt, ocl_resrc->devices, NULL, NULL, &status);

    ocl_resrc->program = ocl_build_from_file(filename, ocl_resrc->context, ocl_resrc->devices, gpu_cnt, &status);
    if (status != CL_SUCCESS || ocl_resrc->program == NULL)
    {
        DBG_PRINT("Failed to build OpenCL program.\n");
        goto RELEASE_CONTEXT;
    }
    else
    {
        /* SUCCESS */
        return status;
    }

/* free resources and return error code */
RELEASE_CONTEXT:
    clReleaseContext(ocl_resrc->context);

FREE_GPUS:
    free(ocl_resrc->devices);
    return status;
}

static cl_int ocl_release(OCLResrc *ocl_resrc)
{
    cl_int status;
    status = clReleaseProgram(ocl_resrc->program);
    status = clReleaseContext(ocl_resrc->context);
    free(ocl_resrc->devices);
    memset(ocl_resrc, 0, sizeof(OCLResrc));
    return status;
}
#endif /* OCL_UTILS_H */
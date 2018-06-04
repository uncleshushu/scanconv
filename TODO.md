# TODO

## Log

1. turn on/off debug prints using build options

1. multi level debug print: info/error/success, with/without filename/line number/function name

## Evaluation

Evaluate the result of interpolation.

## Timing

1. CPU

1. OpenCL

Calculate the acceleration ratio.

## Optimization

Use profiler.

### CPU

1. cache/memory

1. instruction

### OpenCL

1. memory

1. work group

1. eliminate redundant setting-ups

    1. 将调用 kernel 前的所有流程封装成一个函数（是否要在这一步创建 command queue ？），在 `main` 函数中先调用该函数；

    1. 需要用到某种插值算法的 OpenCL 实现时，将对应 kernel 的调用流程封装成函数，调用该函数；

        1. 创建 kernel : `clCreateKernel`

        1. 传入数据:
            
            1. `clCreateBuffer` （用 `CL_MEM_COPY_HOST_PTR` ？）
            
            1. `clSetKernelArg`

        1. 创建 command queue : `clCreateCommandQueueWithProperties`

        1. 设置 kernel 的 global work size 和 local work size

        1. 执行 kernel

        1. 读出结果

        1. 释放 command queue ，内存， kernel

        1. 返回结果

    1. 释放 program, context, devices

## 边界处理

目前只是将边界处偏导都设为 0 。

两种思路

1. 将原图拓展后存为新图，插值时用新图做计算，但是在新图中对应的下标要改

1. 将取原图中某个区域的操作封装为函数，在该函数内部处理边界问题。好像更透明一些。

## 预计算

先把原图中每个区域的插值多项式的系数都计算并保存好。

对于双三次插值，每个函数16个系数，原图有大约 200 * 2000 个像素，大概需要 25 MB 内存。

OpenCL 怎么实现？

写两个 kernel，一个预计算系数，一个求值

1. 在 host code 先后调用这两个 kernel

1. kernel 调 kernel
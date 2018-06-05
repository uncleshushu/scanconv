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

#### 预计算

先把原图中每个区域的插值多项式的系数都计算并保存好。

对于双三次插值，每个函数16个系数，原图有大约 200 * 2000 个像素，大概需要 25 MB 内存。

OpenCL 怎么实现？

写两个 kernel，一个预计算系数，一个求值

1. 在 host code 先后调用这两个 kernel

1. kernel 调 kernel

#### memory

1. 使用 local memory (work group) ，先将一个 work group 需要的数据读入到 local memory 。

    如何确定合适的 group size:

    - `clGetDeviceInfo` with `CL_DEVICE_LOCAL_MEM_SIZE` and `CL_DEVICE_MAX_WORK_ITEM_SIZES`
    
    - `clGetKernelWorkGroupInfo` with `CL_KERNEL_WORK_GROUP_SIZE` and `CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE`
    
    - always make sure `group_size` <= `CL_KERNEL_WORK_GROUP_SIZE` and `local_work_size[i]` <= `CL_DEVICE_MAX_WORK_ITEM_SIZES[i]` and    `group_size` * `local_mem_per_item` <= `CL_DEVICE_LOCAL_MEM_SIZE`
     
    - group sizes of multiples of `CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE` are preferred

1. 对 global memory 访问的优化，如减少 bank-conflict 等。

1. 使用 constant memory

    >常量内存使用__constant标识进行描述，常量内存作为全局地址空间的一部分，在运行时可以分配出出相应的缓存空间，利用常量内存可以提高应用访存效率。使用常量地址空间的方式有两种：
    >
    >1. 可以通过数组创建的方式，之后将数组作为参数传入内核中。内核参数描述上，必须指定__constant为对应指针的标识符。
    >
    >1. 内核端声明常量对象，并使用__constant标识对其进行初始化，其属于编译时常量类型。
    >
    >不同架构下，常量内存对应的位置也不同。AMD GPU上常量数据将会保存在缓存上(保存通用数据)。该缓存比L1缓存的延迟还要低，并且能大大减少GPU内的数据传输。根据这种访存模式，其地址可能就直接保存在指令中，就算是释放功能单元也只有极低的寻址延迟。

    >统一地址空间包括全局、局部和私有地址空间。2.0标准中，并未将常量地址划分到同一地址空间中。虽然，常量空间在逻辑上是全局便令的一部分，不过在某些处理器上(特别是图像处理器)，常量数据会映射到特定的硬件单元上，并不能使用指令对其进行动态指定。
    
    经测试，将 `M_inv[16][16]` 由 private memory 改为 constant memory 后，反而变慢了。

#### kernel 内循环

考虑每个 work item 处理多个像素。

#### 分支消除

1. fmin, fmax, if, ternary operator

#### eliminate redundant setting-ups

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

两种思路

1. 将原图拓展后存为新图，插值时用新图做计算，但是在新图中对应的下标要改

1. 将取原图中某个区域的操作封装为函数，在该函数内部处理边界问题。好像更透明一些。

<del>目前只是将边界处偏导都设为 0 。</del> 目前采用用思路 2 。
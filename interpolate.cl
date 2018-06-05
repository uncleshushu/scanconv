
typedef enum {BORDER_REFLECT} border_t;
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
inline int usi_get_patch(float *usi_pixels, int line_cnt, int spl, 
                            int left, int top, 
                            float *patch, int w, int h, 
                            border_t border_type)
{
    // size limit
    if(w < 0 || h < 0 || w > line_cnt || h > spl)
        return -1;
    // border limit
    if(left < -1 || top < -1 || left + w - 1 > line_cnt || top + h -1 > spl)
        return -1;
    
    if(border_type == BORDER_REFLECT)
    {
        int x, y;
        // int c;    
        for(int i = 0; i < h; ++i)
        {
            y = top + i;
            // c = y < 0;
            // y = c * (-y - 1) + (!c) * y;
            // c = y > spl - 1;
            // y = c * (2*spl - y - 1) + (!c) * y;
            y = y < 0 ? -y-1 : y;
            y = y > spl - 1 ? 2*spl - y - 1 : y;
            for(int j = 0; j < w; ++j)
            {
                x = left + j;
                // c = x < 0;
                // x = c * (-x -1) + (!c) * x;
                // c = x > line_cnt - 1;
                // x = c * (2*line_cnt - x - 1) + (!c) * x;
                x = x < 0 ? -x -1 : x;
                x = x > line_cnt - 1 ? 2*line_cnt - x - 1 : x;
                patch[i*w + j] = usi_pixels[y*line_cnt + x]; 
            }
        } 
    }
    return 0;
}


inline float itp_linear(float y1, float y2, float dx)
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
inline float itp_catmull_rom_spline(float y[4], float dx)
{
    float a[4];
    a[3] = -0.5f * y[0] + 1.5f * y[1] - 1.5f * y[2] + 0.5f * y[3];
    a[2] = y[0] - 2.5f * y[1] + 2.0f * y[2] - 0.5f * y[3]; 
    a[1] = -0.5f * y[0] + 0.5f * y[2];
    a[0] = y[1];
    return fma(fma(fma(a[3], dx, a[2]), dx, a[1]), dx, a[0]);
}



__kernel void nearest(__global float *usi_pixels,    
                            float radius, float angle,  
                            int spl, int line_cnt, 
                            float R, float center2top,
                            float real_w, 
                            float s_interval, float a_interval, 
                            __global unsigned char *gi_pixels)
{
    size_t gi_height = get_global_size(0);
    size_t gi_width = get_global_size(1);
    size_t i = get_global_id(0);
    size_t j = get_global_id(1);

    float real_y = center2top + i*s_interval;
    float real_x = j*s_interval - real_w/2;
    float theta = atan(real_x/real_y);
    float rho = sqrt(real_x * real_x + real_y * real_y);
    if( theta >= -angle && theta <= angle 
        && rho <= R && rho >= radius)
    {
        int j_0 = (int)((theta + angle) / a_interval);
        int i_0 = (int)((rho - radius) / s_interval);
        gi_pixels[i*gi_width + j] = (unsigned char)round(usi_pixels[i_0*line_cnt + j_0]);
    }
}

__kernel void bilinear(__global float *usi_pixels,    
                            float radius, float angle,  
                            int spl, int line_cnt, 
                            float R, float center2top,
                            float real_w, 
                            float s_interval, float a_interval, 
                            __global unsigned char *gi_pixels)
{
    size_t gi_height = get_global_size(0);
    size_t gi_width = get_global_size(1);
    size_t i = get_global_id(0);
    size_t j = get_global_id(1);

    float real_y = center2top + i*s_interval;
    float real_x = j*s_interval - real_w/2;
    float theta = atan(real_x/real_y);
    float rho = sqrt(real_x * real_x + real_y * real_y);
    if( theta >= -angle && theta <= angle 
        && rho <= R && rho >= radius)
    {
        float u = (theta + angle) / a_interval;
        float v = (rho - radius) / s_interval;
        int left = (int)u;
        int right = left + 1;
        int above = (int)v;
        int below = above + 1;

        float pixel_above_left = usi_pixels[above*line_cnt + left];
        float pixel_above_right = usi_pixels[above*line_cnt + right];
        float pixel_below_left = usi_pixels[below*line_cnt + left];
        float pixel_below_right = usi_pixels[below*line_cnt + right];

        float pixel_above = (u - left)*pixel_above_right - (u - right)*pixel_above_left;
        float pixel_below = (u - left)*pixel_below_right - (u - right)*pixel_below_left;  
        
        float pixel_bilinear = (v - above)*pixel_below - (v - below)*pixel_above;
        unsigned char pixel = (unsigned char)round(pixel_bilinear);
        gi_pixels[i*gi_width + j] = pixel;
    }
}

inline float mydot(__constant float a[], const float b[], const size_t n)
{
    float sum = 0;
    for(size_t i = 0; i < n; ++i)
        sum = fma(a[i], b[i], sum);
    return sum;
}
__constant float M_inv[16][16] = 
{
    1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
    -3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0,
    2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
    0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0,
    0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0,
    -3,0,3,0,0,0,0,0,-2,0,-1,0,0,0,0,0,
    0,0,0,0,-3,0,3,0,0,0,0,0,-2,0,-1,0,
    9,-9,-9,9,6,3,-6,-3,6,-6,3,-3,4,2,2,1,
    -6,6,6,-6,-3,-3,3,3,-4,4,-2,2,-2,-2,-1,-1,
    2,0,-2,0,0,0,0,0,1,0,1,0,0,0,0,0,
    0,0,0,0,2,0,-2,0,0,0,0,0,1,0,1,0,
    -6,6,6,-6,-4,-2,4,2,-3,3,-3,3,-2,-1,-2,-1,
    4,-4,-4,4,2,2,-2,-2,2,-2,2,-2,1,1,1,1
};
__kernel void bicubic(__global float *usi_pixels,    
                            float radius, float angle,  
                            int spl, int line_cnt, 
                            float R, float center2top,
                            float real_w, 
                            float s_interval, float a_interval, 
                            __global unsigned char *gi_pixels)
{
    size_t gi_height = get_global_size(0);
    size_t gi_width = get_global_size(1);
    size_t i = get_global_id(0);
    size_t j = get_global_id(1);
    // const float M_inv[16][16] = 
    // {
    //     1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    //     0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
    //     -3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0,
    //     2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0,
    //     0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
    //     0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
    //     0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0,
    //     0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0,
    //     -3,0,3,0,0,0,0,0,-2,0,-1,0,0,0,0,0,
    //     0,0,0,0,-3,0,3,0,0,0,0,0,-2,0,-1,0,
    //     9,-9,-9,9,6,3,-6,-3,6,-6,3,-3,4,2,2,1,
    //     -6,6,6,-6,-3,-3,3,3,-4,4,-2,2,-2,-2,-1,-1,
    //     2,0,-2,0,0,0,0,0,1,0,1,0,0,0,0,0,
    //     0,0,0,0,2,0,-2,0,0,0,0,0,1,0,1,0,
    //     -6,6,6,-6,-4,-2,4,2,-3,3,-3,3,-2,-1,-2,-1,
    //     4,-4,-4,4,2,2,-2,-2,2,-2,2,-2,1,1,1,1
    // };
    float real_y = center2top + i*s_interval;
    float real_x = j*s_interval - real_w/2;
    float theta = atan(real_x/real_y);
    float rho = sqrt(real_x * real_x + real_y * real_y);
    if( theta >= -angle && theta <= angle 
        && rho <= R && rho >= radius)
    {
        float u = (theta + angle) / a_interval;
        float v = (rho - radius) / s_interval;
        int p = (int)u; /* column index */
        int q = (int)v; /* row index */  
        float dx = u - p;
        float dy = v - q; 

        /* derivs = [f(0,0) f(1,0) f(0,1) f(1,1) fx(0,0) ... fy(0,0) ...fxy(0,0)...]^T */
        float derivs[16] = {0.0};
        /* coefficients of the bicubic interpolation polynomial
        * alpha = M_inv * derivs
        * alpha = [a00 a10 a20 a30 a01 a11 a21 a31 ... a33]^T;
        */
        float alpha[16];

        float (*f)[2] = derivs;
        float (*fx)[2] = derivs + 4;
        float (*fy)[2] = derivs + 8;
        float (*fxy)[2] = derivs + 12;
        
        // f(1,0) = f[0][1], f(0,1) = f[1][0]
        f[0][0] = usi_pixels[q*line_cnt + p];
        f[0][1] = usi_pixels[q*line_cnt + (p+1)];
        f[1][0] = usi_pixels[(q+1)*line_cnt + p];
        f[1][1] = usi_pixels[(q+1)*line_cnt + (p+1)];
        
        // if(q > 0 && q < spl - 1 && p > 0 && p < line_cnt -1)
        {
            float f_4x4[4][4];
            // for(int t = 0; t < 4; ++t)
            //     for(int s = 0; s < 4; ++s)
            //         f_4x4[t][s] = usi_pixels[(q-1+t)*line_cnt + (p-1+s)];
            usi_get_patch(usi_pixels, line_cnt, spl, p-1, q-1, f_4x4, 4, 4, BORDER_REFLECT);

            for(int t = 0; t < 2; ++t)
                for(int s = 0; s < 2; ++s)
                {
                    // fx(s,t) = (f(s+1,t) - f(s-1,t)) / 2
                    fx[t][s] = 0.5f * (f_4x4[1+t][1+s+1] - f_4x4[1+t][1+s-1]);
                    // fy(s,t) = (f(s,t+1) - f(s,t-1)) / 2
                    fy[t][s] = 0.5f * (f_4x4[1+t+1][1+s] - f_4x4[1+t-1][1+s]);
                }
            for(int t = 0; t < 2; ++t)
                for(int s = 0; s < 2; ++s)
                    // fxy(s,t) = (fx(s,t+1) - fx(s,t-1)) / 2
                    //          = (f(s+1,t+1) + f(s-1,t-1) - f(s+1,t-1) - f(s-1, t+1)) / 4
                    fy[t][s] = 0.25f * (fx[1+t+1][1+s+1] + f_4x4[1+t-1][1+s-1] 
                                        - f_4x4[1+t-1][1+s+1] - f_4x4[1+t+1][1+s-1]);
        }

        for(size_t n = 0; n < 16; ++n)
            alpha[n] =mydot(M_inv[n], derivs, 16);

        // aij = a[j][i]
        float (*a)[4] = alpha;

        float result = 0.0f;
        for(int n = 3; n >= 0; --n)
        {
            float x_part = 0.0f;    
            for(int m = 3; m >= 0; --m)
                x_part =  fma(x_part, dx, a[n][m]);
            
            result = fma(result, dy, x_part);
        }

        result = fmax(0.0f, result);
        result = fmin(255.0f, result); 
        gi_pixels[i*gi_width + j] = (unsigned char)round(result);
    }
}


__kernel void col_linear_row_cubic(__global float *usi_pixels,    
                            float radius, float angle,  
                            int spl, int line_cnt, 
                            float R, float center2top,
                            float real_w, 
                            float s_interval, float a_interval, 
                            __global unsigned char *gi_pixels)
{
    size_t gi_height = get_global_size(0);
    size_t gi_width = get_global_size(1);
    size_t i = get_global_id(0);
    size_t j = get_global_id(1);

    float real_y = center2top + i*s_interval;
    float real_x = j*s_interval - real_w/2;
    float theta = atan(real_x/real_y);
    float rho = sqrt(real_x * real_x + real_y * real_y);
    if( theta >= -angle && theta <= angle 
        && rho <= R && rho >= radius)
    {
        float u = (theta + angle) / a_interval;
        float v = (rho - radius) / s_interval;
        int p = (int) u;
        int q = (int) v;
        float dx = u - p;
        float dy = v - q;

        float f_2x4[2][4];
        usi_get_patch(usi_pixels, line_cnt, spl, p-1, q-1, f_2x4, 4, 2, BORDER_REFLECT);
        // linear interpolation in y direction
        float lin_y[4];
        for(int n = 0; n < 4; ++n)
            lin_y[n] = itp_linear(f_2x4[0][n], f_2x4[1][n], dy);
        
        // cubic interpolation in x direction
        float cub_x = itp_catmull_rom_spline(lin_y, dx);
        cub_x = cub_x < 0.0f ? 0.0f : cub_x;
        cub_x = cub_x > 255.0f ? 255.0f : cub_x;
        // cub_x = fmax(0.0f, cub_x);
        // cub_x = fmin(255.0f, cub_x); 
        gi_pixels[i*gi_width + j] = (unsigned char) cub_x;
    }
}

__kernel void bicubic_catmull_rom_spline(__global float *usi_pixels,    
                            float radius, float angle,  
                            int spl, int line_cnt, 
                            float R, float center2top,
                            float real_w, 
                            float s_interval, float a_interval, 
                            __global unsigned char *gi_pixels)
{
    size_t gi_height = get_global_size(0);
    size_t gi_width = get_global_size(1);
    size_t i = get_global_id(0);
    size_t j = get_global_id(1);

    float real_y = center2top + i*s_interval;
    float real_x = j*s_interval - real_w/2;
    float theta = atan(real_x/real_y);
    float rho = sqrt(real_x * real_x + real_y * real_y);
    if( theta >= -angle && theta <= angle 
        && rho <= R && rho >= radius)
    {
        float u = (theta + angle) / a_interval;
        float v = (rho - radius) / s_interval;
        int p = (int) u;
        int q = (int) v;
        float dx = u - p;
        float dy = v - q;

        float f_4x4[4][4];
        usi_get_patch(usi_pixels, line_cnt, spl, p-1, q-1, f_4x4, 4, 4, BORDER_REFLECT);
        
        // Catmull-Rom spline interpolation in x direction
        float cub_x[4];
        for(int n = 0; n < 4; ++n)
            cub_x[n] = itp_catmull_rom_spline(f_4x4[n], dx);
        
        // Catmull-Rom spline interpolation in y direction
        float cub_y = itp_catmull_rom_spline(cub_x, dy);
        cub_y = cub_y < 0.0f ? 0.0f : cub_y;
        cub_y = cub_y > 255.0f ? 255.0f : cub_y;
        // cub_y = fmax(0.0f, cub_y);
        // cub_y = fmin(255.0f, cub_y); 
        gi_pixels[i*gi_width + j] = (unsigned char) cub_y;
    }
}
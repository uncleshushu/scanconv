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

inline float mydot(const float a[], const float b[], const size_t n)
{
    float sum = 0;
    for(size_t i = 0; i < n; ++i)
        sum = fma(a[i], b[i], sum);
    return sum;
}
// __constant float M_inv[16][16] = 
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
    const float M_inv[16][16] = 
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
        
        if(q > 0 && q < spl - 1 && p > 0 && p < line_cnt -1)
        {
            float f_4x4[4][4];
            for(int t = 0; t < 4; ++t)
                for(int s = 0; s < 4; ++s)
                    f_4x4[t][s] = usi_pixels[(q-1+t)*line_cnt + (p-1+s)];

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
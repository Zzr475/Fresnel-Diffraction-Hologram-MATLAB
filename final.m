%% Image Read and preprocess
img = imread('input_image3.jpg');% 从文件读取彩色图像
R_img = img(:,:,1); % 提取红色通道
G_img = img(:,:,2); 
B_img = img(:,:,3); 
[M0,N0] = size(R_img); % 获取红色通道图像的尺寸

save_dir = 'Generated_Images'; % 定义保存图片的文件夹名称
if ~exist(save_dir, 'dir') % 检查文件夹是否存在，如果不存在则创建
    mkdir(save_dir);
end
imwrite(img, fullfile(save_dir, 'original_input_image.png'));% 保存原始输入图像

%% Set Basic parameters
N1 = min(M0,N0); % 获取图像宽度和高度的最小值
N = 1024;% 设置全息图采样点数为1024
hologram_R = imresize(R_img,N/4/N1); % 调整红色通道图像大小
hologram_G = imresize(G_img,N/4/N1); 
hologram_B = imresize(B_img,N/4/N1); 
[M1,N1] = size(hologram_R); % 获取调整后红色全息图的尺寸

Xr = zeros(N,N);% 创建N×N零矩阵作为红色通道画布
Xr(N/2-M1/2+1:N/2+M1/2,N/2-N1/2+1:N/2+N1/2) = hologram_R(1:M1,1:N1); % 将红色全息图置于画布中心
hr = 0.630e-3;% 红色波长：0.630毫米
Xg = zeros(N,N); 
Xg(N/2-M1/2+1:N/2+M1/2,N/2-N1/2+1:N/2+N1/2) = hologram_G(1:M1,1:N1);
hg = 0.550e-3; 
Xb = zeros(N,N);
Xb(N/2-M1/2+1:N/2+M1/2,N/2-N1/2+1:N/2+N1/2) = hologram_B(1:M1,1:N1);
hb = 0.470e-3; 
kr = 2*pi/hr; % 计算红色波数
kg = 2*pi/hg;
kb = 2*pi/hb; 
pix = 0.008; % SLM像素宽度：0.008毫米
L = N*pix; % 计算SLM宽度
z0 = 1000; % 全息记录距离：1000毫米
L0r = hr*N*z0/L; % 红色物平面宽度
L0g = hg*N*z0/L;
L0b = hb*N*z0/L; 
L0 = max([L0r,L0g,L0b]); % 取最大波长对应的统一物平面宽度
Yr = double(Xr); % 将红色画布转换为双精度
Yg = double(Xg);
Yb = double(Xb); 
b = rand(N,N)*2*pi;% 生成[0,2π]随机相位矩阵
fr = Yr.*exp(1i.*b);% 为红色通道应用随机相位
fg = Yg.*exp(1i.*b);
fb = Yb.*exp(1i.*b); 
figstrr = strcat('R通道初始物平面宽度=',num2str(L0),'mm');
figstrg = strcat('G通道初始物平面宽度=',num2str(L0),'mm');
figstrb = strcat('B通道初始物平面宽度=',num2str(L0),'mm');
figure(1);
subplot(2,2,1); imshow(img); title('原始图像');
subplot(2,2,2); imshow(Xr,[]); xlabel(figstrr); title('R通道物平面图像');
subplot(2,2,3); imshow(Xg,[]); xlabel(figstrg); title('G通道物平面图像');
subplot(2,2,4); imshow(Xb,[]); xlabel(figstrb); title('B通道物平面图像');
imwrite(mat2gray(Xr), fullfile(save_dir, 'object_plane_R.png'));% 保存 R通道物平面图像
imwrite(mat2gray(Xg), fullfile(save_dir, 'object_plane_G.png'));% 保存 G通道物平面图像
imwrite(mat2gray(Xb), fullfile(save_dir, 'object_plane_B.png'));% 保存 B通道物平面图像

%% Use Fresnel Function
Ufr = Fresnell(fr, N, L0, z0, kr, L, hr); % 计算红色通道物场
Ufg = Fresnell(fg, N, L0, z0, kg, L, hg); 
Ufb = Fresnell(fb, N, L0, z0, kb, L, hb); 
figure(2);
figstrf = strcat('模拟SLM宽度=',num2str(L),'mm');
subplot(2,2,1); imshow(abs(Ufr),[]); xlabel(figstrf); title('R通道SLM平面振幅分布图像');
subplot(2,2,2); imshow(abs(Ufg),[]); xlabel(figstrf); title('G通道SLM平面振幅分布图像');
subplot(2,2,3); imshow(abs(Ufb),[]); xlabel(figstrf); title('B通道SLM平面振幅分布图像');
imwrite(mat2gray(abs(Ufr)), fullfile(save_dir, 'slm_plane_amplitude_R.png'));% 保存 R通道SLM平面振幅分布图像
imwrite(mat2gray(abs(Ufg)), fullfile(save_dir, 'slm_plane_amplitude_G.png'));% 保存 G通道SLM平面振幅分布图像
imwrite(mat2gray(abs(Ufb)), fullfile(save_dir, 'slm_plane_amplitude_B.png'));% 保存 B通道SLM平面振幅分布图像

%% Fresnel Function
function Uf = Fresnell(f, N, L0, z0, k, L, h)
    n = 1:N;% 生成1到N的序列
    x_obj = -L0/2+L0/N*(n-1); y_obj = x_obj; % 计算物平面坐标
    [yy_obj,xx_obj] = meshgrid(y_obj,x_obj);% 生成物平面网格
    x_slm = -L/2+L/N*(n-1);% 计算CCD平面坐标
    y_slm = x_slm; [yy_ccd,xx_ccd] = meshgrid(y_slm,x_slm);% 生成SLM平面网格
    Fresnel = exp(1i*k/2/z0*(xx_obj.^2+yy_obj.^2)); % 计算物平面菲涅耳相位
    f2 = f.*Fresnel; % 应用菲涅耳相位
    Uf = fft2(f2,N,N); % 2D FFT
    Uf = fftshift(Uf);% 将零频移至中心
    phase = exp(1i*k*z0)/(1i*h*z0)*exp(1i*k/2/z0*(xx_ccd.^2+yy_ccd.^2));% 计算CCD平面相位
    Uf = Uf.*phase*(L/N)^2;% 归一化傅里叶变换幅度
end

%% Reference Light Function
function Ur = RL(L0, z0, L, N, Uf, k)
    if k == 2*pi/0.630e-3
        Qx = (1.5)*L0/8/z0;% 红色通道方向余弦
        Qy = Qx;
    elseif k == 2*pi/0.550e-3
        Qx = (1.3095)*L0/8/z0; 
        Qy = Qx;
    else
        Qx = (1.1190)*L0/8/z0; 
        Qy = Qx;
    end
    x_rf = -L/2:L/N:L/2-L/N;
    y_rf = x_rf; % 参考平面坐标
    [X_rf,Y_rf] = meshgrid(x_rf,y_rf);% 生成参考平面网格
    Ar = max(max(abs(Uf))); % 从物场最大值获取参考光幅度
    Ur = Ar*exp(1i*k*(X_rf.*Qx+Y_rf.*Qy));% 计算参考光场
end

%% Interference with Four-Step Phase Shifting
% 四步相位移干涉
Urr_0 = RL(L0, z0, L, N, Ufr, kr); % 红色参考光，0相位
Urr_pi2 = RL(L0, z0, L, N, Ufr, kr) * exp(1i*pi/2); % 红色参考光，π/2相位
Urr_pi = RL(L0, z0, L, N, Ufr, kr) * exp(1i*pi);% 红色参考光，π相位
Urr_3pi2 = RL(L0, z0, L, N, Ufr, kr) * exp(1i*3*pi/2); % 红色参考光，3π/2相位
Urg_0 = RL(L0, z0, L, N, Ufg, kg); 
Urg_pi2 = RL(L0, z0, L, N, Ufg, kg) * exp(1i*pi/2); 
Urg_pi = RL(L0, z0, L, N, Ufg, kg) * exp(1i*pi);
Urg_3pi2 = RL(L0, z0, L, N, Ufg, kg) * exp(1i*3*pi/2); 
Urb_0 = RL(L0, z0, L, N, Ufb, kb); 
Urb_pi2 = RL(L0, z0, L, N, Ufb, kb) * exp(1i*pi/2); 
Urb_pi = RL(L0, z0, L, N, Ufb, kb) * exp(1i*pi);
Urb_3pi2 = RL(L0, z0, L, N, Ufb, kb) * exp(1i*3*pi/2);

% 计算每个相位移的强度
I1r = abs(Ufr + Urr_0).^2; % 红色强度，0相位
I2r = abs(Ufr + Urr_pi2).^2;% 红色强度，π/2相位
I3r = abs(Ufr + Urr_pi).^2;% 红色强度，π相位
I4r = abs(Ufr + Urr_3pi2).^2; % 红色强度，3π/2相位
I1g = abs(Ufg + Urg_0).^2;
I2g = abs(Ufg + Urg_pi2).^2; 
I3g = abs(Ufg + Urg_pi).^2; 
I4g = abs(Ufg + Urg_3pi2).^2; 
I1b = abs(Ufb + Urb_0).^2;
I2b = abs(Ufb + Urb_pi2).^2;
I3b = abs(Ufb + Urb_pi).^2; 
I4b = abs(Ufb + Urb_3pi2).^2;

% 使用四步相位移公式计算复数物波
Uobj_r = ((I1r - I3r) - 1i*(I4r - I2r)) / (4 * max(max(abs(Urr_0))));% 红色物波
Uobj_g = ((I1g - I3g) - 1i*(I4g - I2g)) / (4 * max(max(abs(Urg_0)))); 
Uobj_b = ((I1b - I3b) - 1i*(I4b - I2b)) / (4 * max(max(abs(Urb_0))));

% 将复数物波作为全息图
Ihr = Uobj_r;% 红色复数全息图
Ihg = Uobj_g; 
Ihb = Uobj_b;

% 全息图显示
figure(3)
figstr_inf = strcat('全息图宽度=',num2str(L),'mm');
subplot(2,2,1); imshow(abs(Ihr),[]); xlabel(figstr_inf); title('R通道模拟形成的数字全息图（振幅）');
subplot(2,2,2); imshow(abs(Ihg),[]); xlabel(figstr_inf); title('G通道模拟形成的数字全息图（振幅）');
subplot(2,2,3); imshow(abs(Ihb),[]); xlabel(figstr_inf); title('B通道模拟形成的数字全息图（振幅）');
imwrite(mat2gray(abs(Ihr)), fullfile(save_dir, 'hologram_amplitude_R.png'));% 保存 R通道全息图振幅
imwrite(mat2gray(abs(Ihg)), fullfile(save_dir, 'hologram_amplitude_G.png'));% 保存 G通道全息图振幅
imwrite(mat2gray(abs(Ihb)), fullfile(save_dir, 'hologram_amplitude_B.png'));% 保存 B通道全息图振幅

%% Reconstruction Function
function Uo = Reconstract(Ih, N, L, k, z0, L0, h)
    f0 = double(Ih); 
    In(1:N,1:N) = f0(1:N,1:N); % 将全息图赋值给输入矩阵
    n = 1:N;
    x_re = -L/2+L/N*(n-1); % 计算重建平面坐标
    y_re = x_re; [yy_re,xx_re] = meshgrid(y_re,x_re); % 生成重建平面网格
    Fresnel = exp(-1i*k/2/z0*(xx_re.^2+yy_re.^2)); % 计算重建平面菲涅耳相位（共轭）
    f2 = In.*Fresnel; % 应用菲涅耳相位
    Uf = fft2(f2,N,N); % 执行二维快速傅里叶变换
    Uf = fftshift(Uf)*(L/N)^2; % 将零频移至中心并归一化
    x = -L0/2+L0/N*(n-1); 
    y = x; [yy,xx] = meshgrid(y,x); % Object plane grid
    phase = exp(1i*k*z0)/(1i*h*z0)*exp(1i*k/2/z0*(xx.^2+yy.^2)); 
    Uo = Uf.*phase;% 计算重建物场
    Uo = rot90(Uo, 2); % 将重建场旋转180度以校正方向
end

%% Use Reconstruction Function
Uor = Reconstract(Uobj_r, N, L, kr, z0, L0r, hr); % 重建红色通道
Uog = Reconstract(Uobj_g, N, L, kg, z0, L0g, hg);
Uob = Reconstract(Uobj_b, N, L, kb, z0, L0b, hb); 

%% Show
Gmaxr = max(max(abs(Uor))); % 红色最大振幅
Gminr = min(min(abs(Uor))); % 红色最小振幅
Gmaxg = max(max(abs(Uog)));
Gming = min(min(abs(Uog))); 
Gmaxb = max(max(abs(Uob)));
Gminb = min(min(abs(Uob))); 
Uor_normalized = mat2gray(abs(Uor)); % 归一化红色通道到[0,1]
Uog_normalized = mat2gray(abs(Uog)); 
Uob_normalized = mat2gray(abs(Uob)); 
figstr = strcat('重建物平面宽度=',num2str(L0),'mm');
figure(4); imshow(Uor_normalized); title('R通道重建图像');
figure(5); imshow(Uog_normalized); title('G通道重建图像');
figure(6); imshow(Uob_normalized); title('B通道重建图像');
imwrite(Uor_normalized, fullfile(save_dir, 'reconstructed_image_R.png'));% 保存 R通道重建图像
imwrite(Uog_normalized, fullfile(save_dir, 'reconstructed_image_G.png'));% 保存 G通道重建图像
imwrite(Uob_normalized, fullfile(save_dir, 'reconstructed_image_B.png'));% 保存 B通道重建图像

%% Mix
[height, width] = size(Uor_normalized); % 获取归一化图像尺寸
rgb_image = zeros(height, width, 3); % 创建空的RGB图像
rgb_image(:,:,1) = min(max(Uor_normalized * 255, 0), 255);% 红色通道赋值
rgb_image(:,:,2) = min(max(Uog_normalized * 255, 0), 255); 
rgb_image(:,:,3) = min(max(Uob_normalized * 255, 0), 255);
rgb_image = uint8(rgb_image); % Convert to uint8
figure(7);
imshow(rgb_image);
title('合成的RGB图像');
imwrite(rgb_image, fullfile(save_dir, 'final_mixed_RGB_image.png'));% 保存合成的RGB图像

%% Image Quality Evaluation
[orig_height, orig_width, ~] = size(img);% 获取原始图像尺寸
[rec_height, rec_width, ~] = size(rgb_image);% 获取重建图像原始尺寸（假设 rgb_image 是输入的原始重建图像）
scale_factor_rec = min(orig_height / rec_height, orig_width / rec_width);% 根据重建图像的原始比例缩放
rgb_image_resized = imresize(rgb_image, scale_factor_rec, 'bicubic');

% 裁剪重建图像以聚焦于可见区域
% 根据可见内容调整这些坐标，需确保比例合理
crop_y_start = 23; % 起始行
crop_y_end = 110; % 结束行（根据可见高度调整）
crop_x_start = 8; % 起始列
crop_x_end = 125; % 结束列（根据可见宽度调整）
if (crop_y_end - crop_y_start + 1) / (crop_x_end - crop_x_start + 1) > 1.5 || ...
   (crop_x_end - crop_x_start + 1) / (crop_y_end - crop_y_start + 1) > 1.5
    fprintf('警告：裁剪区域比例失调，请调整 crop_y_end 或 crop_x_end。\n');% 检查裁剪区域比例
end
rgb_image_cropped = rgb_image_resized(crop_y_start:crop_y_end, crop_x_start:crop_x_end, :);% 裁剪重建图像
[crop_height, crop_width, ~] = size(rgb_image_cropped);% 获取裁剪后重建图像的尺寸

% 计算原始图像的缩放比例（基于裁剪后重建图像的尺寸）
diag_orig = sqrt(orig_height^2 + orig_width^2); % 原始图像对角线长度
diag_crop = sqrt(crop_height^2 + crop_width^2); % 裁剪后重建图像对角线长度
scale_factor_img = min(crop_height / orig_height, crop_width / orig_width);% 计算原始图像缩放比例
img_resized = imresize(img, scale_factor_img, 'bicubic');% 调整原始图像大小
[resized_height, resized_width, ~] = size(img_resized);% 获取缩放后图像的尺寸
canvas = zeros(crop_height, crop_width, 3, 'uint8');% 创建一个与裁剪后重建图像相同尺寸的黑色画布

% 计算偏移量，确保图像居中
canvas_offset_y = max(0, floor((crop_height - resized_height) / 2));
canvas_offset_x = max(0, floor((crop_width - resized_width) / 2));

% 将调整大小后的原始图像放置到画布上
if canvas_offset_y + resized_height <= crop_height && canvas_offset_x + resized_width <= crop_width
    canvas(canvas_offset_y + 1:canvas_offset_y + resized_height, ...
           canvas_offset_x + 1:canvas_offset_x + resized_width, :) = img_resized;% 将调整后的原始图像置于画布
else
    fprintf('警告：缩放后图像尺寸 (%d x %d) 超出画布尺寸 (%d x %d)，进行裁剪以适应。\n', ...
            resized_width, resized_height, crop_width, crop_height);
    % 裁剪图像以适应画布
    y_start = max(1, floor((resized_height - crop_height) / 2) + 1);% 计算裁剪起始行
    y_end = min(resized_height, y_start + crop_height - 1);% 计算裁剪起始行
    x_start = max(1, floor((resized_width - crop_width) / 2) + 1);
    x_end = min(resized_width, x_start + crop_width - 1);
    canvas = img_resized(y_start:y_end, x_start:x_end, :);
end

% 计算每个通道的均方误差 (MSE)
mse_r = mean((double(canvas(:,:,1)) - double(rgb_image_cropped(:,:,1))).^2, 'all');% 计算红色通道均方误差
mse_g = mean((double(canvas(:,:,2)) - double(rgb_image_cropped(:,:,2))).^2, 'all');
mse_b = mean((double(canvas(:,:,3)) - double(rgb_image_cropped(:,:,3))).^2, 'all');
mse_avg = (mse_r + mse_g + mse_b) / 3; % 平均 MSE

max_pixel_value = 255; % 对于 uint8 图像
psnr_r = 10 * log10((max_pixel_value^2) / mse_r);% 计算每个通道的峰值信噪比 (PSNR)
psnr_g = 10 * log10((max_pixel_value^2) / mse_g);
psnr_b = 10 * log10((max_pixel_value^2) / mse_b);
psnr_avg = (psnr_r + psnr_g + psnr_b) / 3; % 平均 PSNR

ssim_r = ssim(canvas(:,:,1), rgb_image_cropped(:,:,1));% 计算每个通道的结构相似性指数 (SSIM)
ssim_g = ssim(canvas(:,:,2), rgb_image_cropped(:,:,2));
ssim_b = ssim(canvas(:,:,3), rgb_image_cropped(:,:,3));
ssim_avg = (ssim_r + ssim_g + ssim_b) / 3; % 平均 SSIM

fprintf('像质评价指标 (画布上的调整大小原始图像 vs 裁剪后重建图像):\n');% 显示结果
fprintf('MSE (R/G/B/平均): %.2f / %.2f / %.2f / %.2f\n', mse_r, mse_g, mse_b, mse_avg);
fprintf('PSNR (R/G/B/平均): %.2f dB / %.2f dB / %.2f dB / %.2f dB\n', psnr_r, psnr_g, psnr_b, psnr_avg);
fprintf('SSIM (R/G/B/平均): %.4f / %.4f / %.4f / %.4f\n', ssim_r, ssim_g, ssim_b, ssim_avg);

figure(8);% 显示画布上的原始图像和裁剪后重建图像并排
subplot(1,2,1); imshow(canvas); title('原始图像 (缩放后置于画布)');
subplot(1,2,2); imshow(rgb_image_cropped); title('重建图像 (裁剪后)');


%本代码基于JackHCC的项目进行优化改进，原项目地址：https://github.com/JackHCC/Computer-Generated-Hologram/blob/5ce72697d2ee0c0e4b94a58bde76eb9e428ac283/Matlab/offaxis_interference_hologram.m#L2

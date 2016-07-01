function r = PSNR(x1, x2)

error = x1 - x2;

r = 10 * log10(255^2 / mean(error(:).^2));

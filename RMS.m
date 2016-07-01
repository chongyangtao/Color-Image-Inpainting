function r = RMS(x1, x2)

error = x1 - x2;

r = sqrt(mean(error(:).^2));

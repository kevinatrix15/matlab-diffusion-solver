function [diff] = getDiffusivity(T)

  baseDiff = 1.1234e-4;
  fullDiff = 7e-4;
  rampT = 700;
  slope = (fullDiff - baseDiff)/rampT;

  if (T < rampT)
    diff = slope*T + baseDiff;
  else
    diff = fullDiff;
  end
end

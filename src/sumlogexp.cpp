#include <math.h>

double logsumexp(double nums[], int ct) {
  double max_exp = nums[0], sum = 0.0;
  int i;

  for (i = 1 ; i < ct ; i++)
    if (nums[i] > max_exp)
      max_exp = nums[i];

  for (i = 0; i < ct ; i++)
    sum += exp(nums[i] - max_exp);

  return log(sum) + max_exp;
}

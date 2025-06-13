#include "../ILFT.h"

double ILFT::trapz(double *f, int lower_bound, int upper_bound)
{
    // If the two bounds are the same then return 0
    if (lower_bound == upper_bound)
        return 0;

    // Initialize
    int k_low = lower_bound;
    int k_up = upper_bound;
    double out_val;

    // Add the start and end values
    out_val = 0.5 * f[k_low] + 0.5 * f[k_up];

    // Loop through the rest of the interval
    for (int k1 = k_low + 1; k1 < k_up; ++k1)
    {
        out_val += f[k1];
    }

    // Return output
    return out_val;
}

double ILFT::simps1(double *f, int lower_bound, int upper_bound)
{
    // Adjust for small interval
    if (lower_bound == upper_bound)
        return 0;
    if (lower_bound - upper_bound == 1)
        return (f[lower_bound] + f[upper_bound]) / 2;

    // Declare variables
    int k_low, k_up, k_up0, oe_lower;
    bool odd_intervals;

    // Initialize the output
    double out_val;

    // Check the intervals
    oe_lower = lower_bound % 2;
    odd_intervals = upper_bound % 2 != lower_bound % 2;

    // Limits of integration (Simpson's)
    k_up = upper_bound;
    k_up0 = k_up;
    k_low = lower_bound;

    // If odd number of intervals, then adjust bound
    if (odd_intervals)
        k_up0 = k_up - 1;

    // Add the start and end values
    out_val = f[k_low] + f[k_up0];

    // Loop through the rest of the interval
    for (int k1 = k_low + 1; k1 < k_up0; ++k1)
    {
        out_val += f[k1] * (2.0 + 2.0 * (k1 % 2 != oe_lower));
    }

    // Divide by 3 per Simpsons rule
    out_val /= 3.0;

    // If the final point was cut off to satisfy even intervals, trapezoid the last point
    if (odd_intervals)
    {
        out_val += (f[k_up] + f[k_up0]) * 0.5;
    }

    // Return output
    return out_val;
}
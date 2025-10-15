#include <cmath>
#include <corecrt_math_defines.h>
inline double XYGetDistance(double x1, double y1, double x2, double y2) {
    return sqrt(pow(x1 - x2, 2) +
                pow(y1 - y2, 2));
}
inline double LatLongGetDistance(double lat1, double lon1, double lat2, double lon2) {
    // Haversine formula
    const double r = 6378.388; // Радиус земли
    const double p = M_PI / 180;

    const double a = 0.5 - cos((lat2 - lat1) * p) / 2
                     + cos(lat1 * p) * cos(lat2 * p) *
                       (1 - cos((lon2 - lon1) * p)) / 2;

    return 2 * r * asin(sqrt(a));
}
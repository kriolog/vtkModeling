#ifndef __VectorMath_h__
#define __VectorMath_h__

#include <math.h>
class VectorMath {
public:
    VectorMath(void);
    ~VectorMath(void);
    inline static void CrossProduct(const double p1[3], const double p2[3], double result[3]) {
        result[0] = p1[1] * p2[2] - p1[2] * p2[1];
        result[1] = p1[2] * p2[0] - p1[0] * p2[2];
        result[2] = p1[0] * p2[1] - p1[1] * p2[0];
    }
    inline static double DotProduct(const double p1[3], const double p2[3]) {
        return    p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2];
    }
    inline static double Length(double p[3]) {
        return sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
    }
    inline static double Area(double v[3]) {
        return (v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    }

    inline static void pts2Vec(double* p1, double* p2, double* vec) {
        vec[0] = p1[0]-p2[0];
        vec[1] = p1[1]-p2[1];
        vec[2] = p1[2]-p2[2];
    }
};

#endif

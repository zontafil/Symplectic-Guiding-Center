#ifndef COORD_UTILS_H
#define COORD_UTILS_H

namespace EMFields {
    Vector3d cyl2cart(Vector3d v, Vector3d x) {
        realnum r = sqrt(x[0]*x[0] + x[1]*x[1]);
        Vector3d ret;
        ret[0] = (v[0]*x[0] - v[1]*x[1]) / r;
        ret[1] = (v[0]*x[1] + v[1]*x[0]) / r;
        ret[2] = v[2];
        return ret;
    }
}

#endif
#ifndef EQDSKREADER_H
#define EQDSKREADER_H
#include <string>

using namespace std;

struct eqdsk {
    int nr;
    int nz;
    double rdim;
    double zdim;
    double rcentr;
    double rleft;
    double zmid;
    double rmaxis;
    double zmaxis;
    double simag;
    double sibry;
    double bcentr;
    double current;
    double* fpol;
    double* pres;
    double* ffprim;
    double* pprime;
    double* psirz;
    double* qpsi;
    int nbbbs;
    int limitr;
    double* rbbbs;
    double* zbbbs;
    double* rlim;
    double* zlim;
    double sepmaxz;
    double sepminz;
    double r_min;
    double r_max;
    double r_grid;
    double z_min;
    double z_max;
    double z_grid;
};

eqdsk readEqdskFile(string filename);

#endif
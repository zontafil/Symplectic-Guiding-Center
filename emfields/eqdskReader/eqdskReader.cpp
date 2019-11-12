#include <fstream>
#include <vector>
#include <sstream>
#include <limits>
#include "./eqdskReader.h"
#include <string>

using namespace std;



void fileSkipNValues(istream& file, int n) {
    string s;
    for (int i = 0; i < n; i++) {
        file >> s;
    }
}

vector<string> splitStringSpace(string str) {
    string buf;                 // Have a buffer string
    stringstream ss(str);       // Insert the string into a stream

    vector<string> tokens; // Create vector to hold our words

    while (ss >> buf)
        tokens.push_back(buf);

    return tokens;
}

eqdsk readEqdskFile(string filename) {
    double inf = numeric_limits<double>::infinity();
    eqdsk ret;

    // open the file
    ifstream file;
    file.open(filename);

    // get nr and nz form the header
    string header;
    getline(file, header);
    vector<string> headerArray = splitStringSpace(header);
    int headerSize = headerArray.size();
    stringstream nr(headerArray[headerSize - 2]);
    nr >> ret.nr;
    stringstream nz(headerArray[headerSize - 1]);
    nz >> ret.nz;

    // read (neqdsk,2020) rdim,zdim,rcentr,rleft,zmid
    file >> ret.rdim;
    file >> ret.zdim;
    file >> ret.rcentr;
    file >> ret.rleft;
    file >> ret.zmid;

    // read (neqdsk,2020) rmaxis,zmaxis,simag,sibry,bcentr
    file >> ret.rmaxis;
    file >> ret.zmaxis;
    file >> ret.simag;
    file >> ret.sibry;
    file >> ret.bcentr;
    file >> ret.current;
    fileSkipNValues(file, 9);
    ret.fpol = (double*)malloc(ret.nr*sizeof(double));
    for (int i=0; i < ret.nr; i++) {
        file >> ret.fpol[i];
    }
    ret.pres = (double*)malloc(ret.nr*sizeof(double));
    for (int i=0; i < ret.nr; i++) {
        file >> ret.pres[i];
    }
    ret.ffprim = (double*)malloc(ret.nr*sizeof(double));
    for (int i=0; i < ret.nr; i++) {
        file >> ret.ffprim[i];
    }
    ret.pprime = (double*)malloc(ret.nr*sizeof(double));
    for (int i=0; i < ret.nr; i++) {
        file >> ret.pprime[i];
    }
    ret.psirz = (double*)malloc(ret.nr*ret.nz*sizeof(double));
    for (int i=0; i < ret.nr * ret.nz; i++) {
        file >> ret.psirz[i];
    }
    ret.qpsi = (double*)malloc(ret.nr*sizeof(double));
    for (int i=0; i < ret.nr; i++) {
        file >> ret.qpsi[i];
    }
    file >> ret.nbbbs;
    file >> ret.limitr;

    ret.rbbbs = (double*)malloc(ret.nbbbs*sizeof(double));
    ret.zbbbs = (double*)malloc(ret.nbbbs*sizeof(double));
    for (int i=0; i < 2*ret.nbbbs; i++) {
        if ((i % 2) == 0) {
            file >> ret.rbbbs[i/2];
        } else {
            file >> ret.zbbbs[(i-1)/2];
        }
    }

    ret.rlim = (double*)malloc(ret.limitr*sizeof(double));
    ret.zlim = (double*)malloc(ret.limitr*sizeof(double));
    for (int i=0; i < 2*ret.limitr; i++) {
        if ((i % 2) == 0) {
            file >> ret.rlim[i/2];
        } else {
            file >> ret.zlim[(i-1)/2];
        }
    }

    // Set up the
    // minimum and maximum in z for the separatrix curve
    // for the quick check to evaluate the fpol function
    ret.sepmaxz = -inf;
    ret.sepminz = inf;
    for (int i=0; i < ret.nbbbs; i++) {
        if (ret.zbbbs[i] > ret.sepmaxz) {
            ret.sepmaxz = ret.zbbbs[i];
        }
        if (ret.zbbbs[i] < ret.sepminz) {
            ret.sepminz = ret.zbbbs[i];
        }
    }

    ret.r_min = ret.rleft;
    ret.r_max = ret.r_min + ret.rdim;
    ret.r_grid = (ret.r_max - ret.r_min) / (ret.nr - 1);
    ret.z_min = ret.zmid - ret.zdim / 2.0;
    ret.z_max = ret.zmid + ret.zdim / 2.0;
    ret.z_grid = (ret.z_max - ret.z_min) / (ret.nz - 1);

    return ret;
}
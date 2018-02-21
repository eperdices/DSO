/**
 *
 *  Copyright (C) 2018 Eduardo Perdices <eperdices at gsyc dot es>
 *
 *  The following code is a derivative work of the code from the dso project,
 *  which is licensed under the GNU Public License, version 3. This code therefore
 *  is also licensed under the terms of the GNU Public License, version 3.
 *  For more information see <http://vision.in.tum.de/dso>.
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <sstream>
#include <fstream>
#include <iostream>
#include <iterator>
#include "util/Undistort.h"
#include "util/settings.h"
#include "util/globalFuncs.h"
#include "IOWrapper/ImageDisplay.h"
#include "IOWrapper/ImageRW.h"

namespace dso {

PhotometricUndistorter::PhotometricUndistorter(std::string file, std::string noiseImage,
                                               std::string vignetteImage, int w_, int h_) {
    valid=false;
    vignetteMap=0;
    vignetteMapInv=0;
    w = w_;
    h = h_;
    output = new ImageAndExposure(w,h);
    if(file=="" || vignetteImage=="") {
        printf("NO PHOTOMETRIC Calibration!\n");
    }

    // read G.
    std::ifstream f(file.c_str());
    printf("Reading Photometric Calibration from file %s\n",file.c_str());
    if (!f.good()) {
        printf("PhotometricUndistorter: Could not open file!\n");
        return;
    }

    {
        std::string line;
        std::getline( f, line );
        std::istringstream l1i( line );
        std::vector<float> Gvec = std::vector<float>(std::istream_iterator<float>(l1i), std::istream_iterator<float>());

        GDepth = Gvec.size();

        if(GDepth < 256) {
            printf("PhotometricUndistorter: invalid format! got %d entries in first line, expected at least 256!\n",(int)Gvec.size());
            return;
        }

        for(int i=0;i<GDepth;i++)
            G[i] = Gvec[i];

        for(int i=0;i<GDepth-1;i++) {
            if(G[i+1] <= G[i]) {
                printf("PhotometricUndistorter: G invalid! it has to be strictly increasing, but it isnt!\n");
                return;
            }
        }

        float min=G[0];
        float max=G[GDepth-1];
        for(int i=0;i<GDepth;i++)
            G[i] = 255.0 * (G[i] - min) / (max-min);            // make it to 0..255 => 0..255.
    }

    if(setting_photometricCalibration==0) {
        for(int i=0;i<GDepth;i++)
            G[i] = 255.0f*i/(float)(GDepth-1);
    }

    printf("Reading Vignette Image from %s\n",vignetteImage.c_str());
    MinimalImage<unsigned short>* vm16 = IOWrap::readImageBW_16U(vignetteImage.c_str());
    MinimalImageB* vm8 = IOWrap::readImageBW_8U(vignetteImage.c_str());
    vignetteMap = new float[w*h];
    vignetteMapInv = new float[w*h];

    if(vm16 != 0) {
        if(vm16->w != w ||vm16->h != h) {
            printf("PhotometricUndistorter: Invalid vignette image size! got %d x %d, expected %d x %d\n",
                    vm16->w, vm16->h, w, h);
            if(vm16!=0) delete vm16;
            if(vm8!=0) delete vm8;
            return;
        }

        float maxV=0;
        for(int i=0;i<w*h;i++)
            if(vm16->at(i) > maxV) maxV = vm16->at(i);

        for(int i=0;i<w*h;i++)
            vignetteMap[i] = vm16->at(i) / maxV;
    } else if(vm8 != 0) {
        if(vm8->w != w ||vm8->h != h) {
            printf("PhotometricUndistorter: Invalid vignette image size! got %d x %d, expected %d x %d\n",
                    vm8->w, vm8->h, w, h);
            if(vm16!=0) delete vm16;
            if(vm8!=0) delete vm8;
            return;
        }

        float maxV=0;
        for(int i=0;i<w*h;i++)
            if(vm8->at(i) > maxV) maxV = vm8->at(i);

        for(int i=0;i<w*h;i++)
            vignetteMap[i] = vm8->at(i) / maxV;
    } else {
        printf("PhotometricUndistorter: Invalid vignette image\n");
        if(vm16!=0) delete vm16;
        if(vm8!=0) delete vm8;
        return;
    }

    if(vm16!=0) delete vm16;
    if(vm8!=0) delete vm8;

    for(int i=0;i<w*h;i++)
        vignetteMapInv[i] = 1.0f / vignetteMap[i];

    printf("Successfully read photometric calibration!\n");
    valid = true;
}
PhotometricUndistorter::~PhotometricUndistorter() {
    if(vignetteMap != 0) delete[] vignetteMap;
    if(vignetteMapInv != 0) delete[] vignetteMapInv;
    delete output;
}

void PhotometricUndistorter::unMapFloatImage(float* image) {
    int wh=w*h;
    for(int i=0;i<wh;i++) {
        float BinvC;
        float color = image[i];

        if(color < 1e-3) {
            BinvC=0.0f;
        } else if(color > GDepth-1.01f) {
            BinvC=GDepth-1.1;
        } else {
            int c = color;
            float a = color-c;
            BinvC=G[c]*(1-a) + G[c+1]*a;
        }

        float val = BinvC;
        if(val < 0) val = 0;
        image[i] = val;
    }
}

template<typename T>
void PhotometricUndistorter::processFrame(T* image_in, float exposure_time, float factor) {
    int wh=w*h;
    float* data = output->image;
    assert(output->w == w && output->h == h);
    assert(data != 0);

    if(!valid || exposure_time <= 0 || setting_photometricCalibration==0) { // disable full photometric calibration.
        for(int i=0; i<wh;i++) {
            data[i] = factor*image_in[i];
        }
        output->exposure_time = exposure_time;
        output->timestamp = 0;
    }
    else {
        for(int i=0; i<wh;i++) {
            data[i] = G[image_in[i]];
        }

        if(setting_photometricCalibration==2) {
            for(int i=0; i<wh;i++)
                data[i] *= vignetteMapInv[i];
        }

        output->exposure_time = exposure_time;
        output->timestamp = 0;
    }

    if(!setting_useExposure)
        output->exposure_time = 1;

}
template void PhotometricUndistorter::processFrame<unsigned char>(unsigned char* image_in, float exposure_time, float factor);
template void PhotometricUndistorter::processFrame<unsigned short>(unsigned short* image_in, float exposure_time, float factor);

Undistort::Undistort(const char* configFileName, bool noprefix) {
    printf("Creating undistorter\n");

    if(noprefix)
        readFromFile(configFileName, 9);
    else
        readFromFile(configFileName, 9, "RadTan ");
}

Undistort::~Undistort() {
    if(remapX != 0) delete[] remapX;
    if(remapY != 0) delete[] remapY;
}

Undistort* Undistort::getUndistorterForFile(std::string configFilename, std::string gammaFilename, std::string vignetteFilename) {
    printf("Reading Calibration from file %s",configFilename.c_str());

    std::ifstream f(configFilename.c_str());
    if (!f.good()) {
        f.close();
        printf(" ... not found. Cannot operate without calibration, shutting down.\n");
        f.close();
        return 0;
    }

    printf(" ... found!\n");
    std::string l1;
    std::getline(f,l1);
    f.close();

    float ic[10];

    Undistort* u;

    // Use RadTan model for 9 parameters.
    if(std::sscanf(l1.c_str(), "%f %f %f %f %f %f %f %f %f",
            &ic[0], &ic[1], &ic[2], &ic[3],
            &ic[4], &ic[5], &ic[6], &ic[7], &ic[8]) == 9) {
        printf("found RadTan (OpenCV) camera model, building rectifier.\n");
        u = new Undistort(configFilename.c_str(), true);
            if(!u->isValid()) { delete u; return 0; }
    } else {
        printf("could not read calib file! exit.");
        exit(1);
    }

    u->loadPhotometricCalibration(gammaFilename, "", vignetteFilename);

    return u;
}

void Undistort::loadPhotometricCalibration(std::string file, std::string noiseImage, std::string vignetteImage) {
    photometricUndist = new PhotometricUndistorter(file, noiseImage, vignetteImage,getOriginalSize()[0], getOriginalSize()[1]);
}

template<typename T>
ImageAndExposure* Undistort::undistort(const MinimalImage<T>* image_raw, float exposure, double timestamp, float factor) const {
    if(image_raw->w != wOrg || image_raw->h != hOrg) {
        printf("Undistort::undistort: wrong image size (%d %d instead of %d %d) \n", image_raw->w, image_raw->h, w, h);
        exit(1);
    }

    photometricUndist->processFrame<T>(image_raw->data, exposure, factor);
    ImageAndExposure* result = new ImageAndExposure(w, h, timestamp);
    photometricUndist->output->copyMetaTo(*result);

    if (!passthrough) {
        float* out_data = result->image;
        float* in_data = photometricUndist->output->image;

        float* noiseMapX=0;
        float* noiseMapY=0;
        if(benchmark_varNoise>0) {
            int numnoise=(benchmark_noiseGridsize+8)*(benchmark_noiseGridsize+8);
            noiseMapX=new float[numnoise];
            noiseMapY=new float[numnoise];
            memset(noiseMapX,0,sizeof(float)*numnoise);
            memset(noiseMapY,0,sizeof(float)*numnoise);

            for(int i=0;i<numnoise;i++) {
                noiseMapX[i] =  2*benchmark_varNoise * (rand()/(float)RAND_MAX - 0.5f);
                noiseMapY[i] =  2*benchmark_varNoise * (rand()/(float)RAND_MAX - 0.5f);
            }
        }

        for(int idx = w*h-1;idx>=0;idx--) {
            // get interp. values
            float xx = remapX[idx];
            float yy = remapY[idx];

            if(benchmark_varNoise>0) {
                float deltax = getInterpolatedElement11BiCub(noiseMapX, 4+(xx/(float)wOrg)*benchmark_noiseGridsize, 4+(yy/(float)hOrg)*benchmark_noiseGridsize, benchmark_noiseGridsize+8 );
                float deltay = getInterpolatedElement11BiCub(noiseMapY, 4+(xx/(float)wOrg)*benchmark_noiseGridsize, 4+(yy/(float)hOrg)*benchmark_noiseGridsize, benchmark_noiseGridsize+8 );
                float x = idx%w + deltax;
                float y = idx/w + deltay;
                if(x < 0.01) x = 0.01;
                if(y < 0.01) y = 0.01;
                if(x > w-1.01) x = w-1.01;
                if(y > h-1.01) y = h-1.01;

                xx = getInterpolatedElement(remapX, x, y, w);
                yy = getInterpolatedElement(remapY, x, y, w);
            }

            if(xx<0) {
                out_data[idx] = 0;
            } else {
                // get integer and rational parts
                int xxi = xx;
                int yyi = yy;
                xx -= xxi;
                yy -= yyi;
                float xxyy = xx*yy;

                // get array base pointer
                const float* src = in_data + xxi + yyi * wOrg;

                // interpolate (bilinear)
                out_data[idx] =  xxyy * src[1+wOrg] + (yy-xxyy) * src[wOrg] + (xx-xxyy) * src[1] + (1-xx-yy+xxyy) * src[0];
            }
        }

        if(benchmark_varNoise>0) {
            delete[] noiseMapX;
            delete[] noiseMapY;
        }

    } else {
        memcpy(result->image, photometricUndist->output->image, sizeof(float)*w*h);
    }

    applyBlurNoise(result->image);

    return result;
}

template ImageAndExposure* Undistort::undistort<unsigned char>(const MinimalImage<unsigned char>* image_raw, float exposure, double timestamp, float factor) const;
template ImageAndExposure* Undistort::undistort<unsigned short>(const MinimalImage<unsigned short>* image_raw, float exposure, double timestamp, float factor) const;

void Undistort::applyBlurNoise(float* img) const {
    if(benchmark_varBlurNoise==0) return;

    int numnoise=(benchmark_noiseGridsize+8)*(benchmark_noiseGridsize+8);
    float* noiseMapX=new float[numnoise];
    float* noiseMapY=new float[numnoise];
    float* blutTmp=new float[w*h];

    if(benchmark_varBlurNoise>0) {
        for(int i=0;i<numnoise;i++) {
            noiseMapX[i] =  benchmark_varBlurNoise  * (rand()/(float)RAND_MAX);
            noiseMapY[i] =  benchmark_varBlurNoise  * (rand()/(float)RAND_MAX);
        }
    }

    float gaussMap[1000];
    for(int i=0;i<1000;i++)
        gaussMap[i] = expf((float)(-i*i/(100.0*100.0)));

    // x-blur.
    for(int y=0;y<h;y++) {
        for(int x=0;x<w;x++) {
            float xBlur = getInterpolatedElement11BiCub(noiseMapX, 4+(x/(float)w)*benchmark_noiseGridsize,
                    4+(y/(float)h)*benchmark_noiseGridsize, benchmark_noiseGridsize+8 );

            if(xBlur < 0.01) xBlur=0.01;

            int kernelSize = 1 + (int)(1.0f+xBlur*1.5);
            float sumW=0;
            float sumCW=0;
            for(int dx=0; dx <= kernelSize; dx++) {
                int gmid = 100.0f*dx/xBlur + 0.5f;
                if(gmid > 900 ) gmid = 900;
                float gw = gaussMap[gmid];

                if(x+dx>0 && x+dx<w) {
                    sumW += gw;
                    sumCW += gw * img[x+dx+y*this->w];
                }

                if(x-dx>0 && x-dx<w && dx!=0) {
                    sumW += gw;
                    sumCW += gw * img[x-dx+y*this->w];
                }
            }

            blutTmp[x+y*this->w] = sumCW / sumW;
        }
    }

    // y-blur.
    for(int x=0;x<w;x++) {
        for(int y=0;y<h;y++) {
            float yBlur = getInterpolatedElement11BiCub(noiseMapY,  4+(x/(float)w)*benchmark_noiseGridsize,
                    4+(y/(float)h)*benchmark_noiseGridsize, benchmark_noiseGridsize+8 );

            if(yBlur < 0.01) yBlur=0.01;

            int kernelSize = 1 + (int)(0.9f+yBlur*2.5);
            float sumW=0;
            float sumCW=0;
            for(int dy=0; dy <= kernelSize; dy++) {
                int gmid = 100.0f*dy/yBlur + 0.5f;
                if(gmid > 900 ) gmid = 900;
                float gw = gaussMap[gmid];

                if(y+dy>0 && y+dy<h) {
                    sumW += gw;
                    sumCW += gw * blutTmp[x+(y+dy)*this->w];
                }

                if(y-dy>0 && y-dy<h && dy!=0) {
                    sumW += gw;
                    sumCW += gw * blutTmp[x+(y-dy)*this->w];
                }
            }
            img[x+y*this->w] = sumCW / sumW;
        }
    }

    delete[] noiseMapX;
    delete[] noiseMapY;
}

void Undistort::readFromFile(const char* configFileName, int nPars, std::string prefix) {
    photometricUndist=0;
    valid = false;
    passthrough=false;
    remapX = 0;
    remapY = 0;
    
    parsOrg = VecX(nPars);

    // read parameters
    std::ifstream infile(configFileName);
    assert(infile.good());

    std::string l1,l2;
    std::getline(infile,l1);
    std::getline(infile,l2);

    // l1 & l2
    if(nPars == 9) { // Radtan model
        char buf[1000];
        snprintf(buf, 1000, "%s%%lf %%lf %%lf %%lf %%lf %%lf %%lf %%lf %%lf %%lf %%lf", prefix.c_str());

        if(std::sscanf(l1.c_str(), buf,
                &parsOrg[0], &parsOrg[1], &parsOrg[2], &parsOrg[3], &parsOrg[4],
                &parsOrg[5], &parsOrg[6], &parsOrg[7], &parsOrg[8]) == 9 &&
                std::sscanf(l2.c_str(), "%d %d", &wOrg, &hOrg) == 2) {
            printf("Input resolution: %d %d\n",wOrg, hOrg);
            printf("In: %s%f %f %f %f %f %f %f %f %f\n",
                    prefix.c_str(),
                    parsOrg[0], parsOrg[1], parsOrg[2], parsOrg[3], parsOrg[4], parsOrg[5], parsOrg[6], parsOrg[7], parsOrg[8]);
        } else {
            printf("Failed to read camera calibration (invalid format?)\nCalibration file: %s\n", configFileName);
            infile.close();
            return;
        }
    } else {
        printf("called with invalid number of parameters.... forgot to implement me?\n");
        infile.close();
        return;
    }

    if(parsOrg[2] < 1 && parsOrg[3] < 1) {
        printf("\n\nFound fx=%f, fy=%f, cx=%f, cy=%f.\n I'm assuming this is the \"relative\" calibration file format,"
             "and will rescale this by image width / height to fx=%f, fy=%f, cx=%f, cy=%f.\n\n",
             parsOrg[0], parsOrg[1], parsOrg[2], parsOrg[3],
             parsOrg[0] * wOrg, parsOrg[1] * hOrg, parsOrg[2] * wOrg - 0.5, parsOrg[3] * hOrg - 0.5 );

        // rescale and substract 0.5 offset.
        // the 0.5 is because I'm assuming the calibration is given such that the pixel at (0,0)
        // contains the integral over intensity over [0,0]-[1,1], whereas I assume the pixel (0,0)
        // to contain a sample of the intensity ot [0,0], which is best approximated by the integral over
        // [-0.5,-0.5]-[0.5,0.5]. Thus, the shift by -0.5.
        parsOrg[0] = parsOrg[0] * wOrg;
        parsOrg[1] = parsOrg[1] * hOrg;
        parsOrg[2] = parsOrg[2] * wOrg - 0.5;
        parsOrg[3] = parsOrg[3] * hOrg - 0.5;
    }

    infile.close();

    w = wOrg;
    h = hOrg;

    remapX = new float[w*h];
    remapY = new float[w*h];

    K.setIdentity();
    K(0,0) = parsOrg[0];
    K(1,1) = parsOrg[1];
    K(0,2) = parsOrg[2];
    K(1,2) = parsOrg[3];

    std::cout << "Camera Matrix:" << std::endl << K << std::endl << std::endl;

    for(int y=0;y<h;y++) {
        for(int x=0;x<w;x++) {
            remapX[x+y*w] = x;
            remapY[x+y*w] = y;
        }
    }

    distortCoordinates(remapX, remapY, remapX, remapY, h*w);

    for(int y=0;y<h;y++) {
        for(int x=0;x<w;x++) {
            // make rounding resistant.
            float ix = remapX[x+y*w];
            float iy = remapY[x+y*w];

            if(ix == 0) ix = 0.001;
            if(iy == 0) iy = 0.001;
            if(ix == wOrg-1) ix = wOrg-1.001;
            if(iy == hOrg-1) ix = hOrg-1.001;

            if(ix > 0 && iy > 0 && ix < wOrg-1 &&  iy < wOrg-1) {
                remapX[x+y*w] = ix;
                remapY[x+y*w] = iy;
            } else {
                remapX[x+y*w] = -1;
                remapY[x+y*w] = -1;
            }
        }
    }

    valid = true;
}

void Undistort::distortCoordinates(float* in_x, float* in_y, float* out_x, float* out_y, int n) const {
    // Radial-Tangential distortion
    float fx = parsOrg[0];
    float fy = parsOrg[1];
    float cx = parsOrg[2];
    float cy = parsOrg[3];
    float k1 = parsOrg[4];
    float k2 = parsOrg[5];
    float p1 = parsOrg[6];
    float p2 = parsOrg[7];
    float k3 = parsOrg[8];

    float ofx = K(0,0);
    float ofy = K(1,1);
    float ocx = K(0,2);
    float ocy = K(1,2);

    for(int i=0; i<n; i++) {
        float x = in_x[i];
        float y = in_y[i];

        float ix = (x - ocx) / ofx;
        float iy = (y - ocy) / ofy;

        // Radial
        float rho2 = ix*ix + iy*iy;
        float rad_dist = 1.0 + k1*rho2 + k2*rho2*rho2 + k3*rho2*rho2*rho2;
        float x_dist = ix*rad_dist;
        float y_dist = iy*rad_dist;

        // Tangential
        x_dist = x_dist + 2*p1*ix*iy + p2*(rho2 + 2*ix*ix);
        y_dist = y_dist + 2*p2*ix*iy + p1*(rho2 + 2*iy*iy);

        out_x[i] = fx*x_dist+cx;
        out_y[i] = fy*y_dist+cy;
    }   
}

}

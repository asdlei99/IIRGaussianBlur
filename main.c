
#include "browse.h"

#define USE_SHELL_OPEN

#define STB_IMAGE_STATIC
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "stb_image.h"
#include "stb_image_write.h"

#include <math.h>
#include <stdio.h>
#include "timing.h"
#include <stdint.h>

#ifndef MIN
#define MIN(a, b)    ( (a) > (b) ? (b) : (a) )
#endif
#ifndef MAX
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#endif


unsigned char *loadImage(const char *filename, int *Width, int *Height, int *Channels) {
    return (stbi_load(filename, Width, Height, Channels, 0));
}


void saveImage(const char *filename, int Width, int Height, int Channels, unsigned char *Output) {

    if (!stbi_write_jpg(filename, Width, Height, Channels, Output, 100)) {
        fprintf(stderr, "save JPEG fail.\n");
        return;
    }
#ifdef USE_SHELL_OPEN
    browse(filename);
#endif
}


void splitpath(const char *path, char *drv, char *dir, char *name, char *ext) {
    const char *end;
    const char *p;
    const char *s;
    if (path[0] && path[1] == ':') {
        if (drv) {
            *drv++ = *path++;
            *drv++ = *path++;
            *drv = '\0';
        }
    } else if (drv)
        *drv = '\0';
    for (end = path; *end && *end != ':';)
        end++;
    for (p = end; p > path && *--p != '\\' && *p != '/';)
        if (*p == '.') {
            end = p;
            break;
        }
    if (ext)
        for (s = end; (*ext = *s++);)
            ext++;
    for (p = end; p > path;)
        if (*--p == '\\' || *p == '/') {
            p++;
            break;
        }
    if (name) {
        for (s = p; s < end;)
            *name++ = *s++;
        *name = '\0';
    }
    if (dir) {
        for (s = path; s < p;)
            *dir++ = *s++;
        *dir = '\0';
    }
}


void CalGaussianCoeff(float sigma, float *a0, float *a1, float *a2, float *a3,
                      float *b1, float *b2, float *cprev, float *cnext) {
    if (sigma < 0.5f)
        sigma = 0.5f;
    float alpha = expf((0.726f) * (0.726f)) / sigma;
    float lamma = expf(-alpha);
    *b2 = expf(-2 * alpha);
    float k = (1 - lamma) * (1 - lamma) / (1 + 2 * alpha * lamma - (*b2));
    *a0 = k;
    *a1 = k * (alpha - 1) * lamma;
    *a2 = k * (alpha + 1) * lamma;
    *a3 = -k * (*b2);
    *b1 = -2 * lamma;
    *cprev = (*a0 + *a1) / (1 + *b1 + *b2);
    *cnext = (*a2 + *a3) / (1 + *b1 + *b2);
}

void gaussianHorizontal(float *bufferPerLine, const float *lpRow,
                        float *lpCol, int width, int Channels,
                        int HeightStep, float a0a1, float a2a3, float b1b2,
                        float cprev, float cnext) {
    float *prev = (float *) calloc(Channels, sizeof(float));
    if (prev == NULL) return;
    int lastXPos = width - 1;
    for (int c = 0; c < Channels; ++c) {
        prev[c] = (lpRow[c] * cprev);
    }
    for (int x = 0; x < width; ++x) {
        for (int c = 0; c < Channels; ++c) {
            prev[c] = ((lpRow[c] * (a0a1)) - (prev[c] * (b1b2)));
            bufferPerLine[c] = prev[c];
        }
        bufferPerLine += Channels;
        lpRow += Channels;
    }
    lpRow -= Channels;
    lpCol += HeightStep * lastXPos;
    bufferPerLine -= Channels;
    for (int c = 0; c < Channels; ++c) {
        prev[c] = (lpRow[c] * cnext);
    }
    for (int x = lastXPos; x >= 0; --x) {
        for (int c = 0; c < Channels; ++c) {
            prev[c] = ((lpRow[c] * (a2a3)) - (prev[c] * (b1b2)));
            bufferPerLine[c] += prev[c];
            lpCol[c] = bufferPerLine[c];
        }
        lpRow -= Channels;
        lpCol -= HeightStep;
        bufferPerLine -= Channels;
    }
    free(prev);
}

void
gaussianVertical(float *bufferPerLine, const float *lpRow, float *lpCol, int height, int Channels, int stride,
                 float a0a1, float a2a3, float b1b2, float cprev,
                 float cnext) {
    float *prev = (float *) calloc(Channels, sizeof(float));
    if (prev == NULL) return;
    int lastYPos = height - 1;
    for (int c = 0; c < Channels; ++c) {
        prev[c] = (lpRow[c] * cprev);
    }
    for (int y = 0; y < height; y++) {
        for (int c = 0; c < Channels; ++c) {
            prev[c] = ((lpRow[c] * a0a1) - (prev[c] * b1b2));
            bufferPerLine[c] = prev[c];
        }
        bufferPerLine += Channels;
        lpRow += Channels;
    }
    lpRow -= Channels;
    bufferPerLine -= Channels;
    lpCol += stride * lastYPos;
    for (int c = 0; c < Channels; ++c) {
        prev[c] = (lpRow[c] * cnext);
    }
    for (int y = lastYPos; y >= 0; y--) {
        for (int c = 0; c < Channels; ++c) {
            prev[c] = ((lpRow[c] * a2a3) - (prev[c] * b1b2));
            bufferPerLine[c] += prev[c];
            lpCol[c] = bufferPerLine[c];
        }
        lpRow -= Channels;
        lpCol -= stride;
        bufferPerLine -= Channels;
    }
    free(prev);
}

void GaussianBlur(float *Input, float *Output, int Width, int Height, int Channels, int Stride, float Sigma) {

    float a0, a1, a2, a3, b1, b2, cprev, cnext;
    CalGaussianCoeff(Sigma, &a0, &a1, &a2, &a3, &b1, &b2, &cprev, &cnext);
    float a0a1 = (a0 + a1);
    float a2a3 = (a2 + a3);
    float b1b2 = (b1 + b2);
    int bufferSizePerThread = (Width > Height ? Width : Height) * Channels;
    float *bufferPerLine = (float *) calloc((bufferSizePerThread + Height * Stride), sizeof(float));
    float *cacheData = bufferPerLine + bufferSizePerThread;
    if (bufferPerLine) {
        for (int y = 0; y < Height; ++y) {
            float *lpRow = Input + Stride * y;
            float *lpCol = cacheData + y * Channels;
            gaussianHorizontal(bufferPerLine, lpRow, lpCol, Width, Channels,
                               Channels * Height, a0a1, a2a3, b1b2, cprev, cnext);
        }
        int HeightStep = Height * Channels;
        for (int x = 0; x < Width; ++x) {
            float *lpCol = Output + x * Channels;
            float *lpRow = cacheData + HeightStep * x;
            gaussianVertical(bufferPerLine, lpRow, lpCol, Height, Channels, Stride,
                             a0a1, a2a3, b1b2, cprev, cnext);
        }
        free(bufferPerLine);
    }
}


int main(int argc, char **argv) {
    printf("IIR Gaussian Blur Filter Implementation In C\n ");
    printf("blog:http://cpuimage.cnblogs.com/ \n ");
    if (argc < 2) {
        printf("usage: %s   image \n ", argv[0]);
        printf("eg: %s   d:\\image.jpg \n ", argv[0]);
        return (0);
    }
    char *in_file = argv[1];
    char drive[3];
    char dir[256];
    char fname[256];
    char ext[256];
    char out_file[1024];
    splitpath(in_file, drive, dir, fname, ext);
    sprintf(out_file, "%s%s%s_out.jpg", drive, dir, fname);

    int Width = 0;
    int Height = 0;
    int Channels = 0;
    unsigned char *inputImage = NULL;
    double startTime = now();
    inputImage = loadImage(in_file, &Width, &Height, &Channels);
    double nLoadTime = calcElapsed(startTime, now());
    printf("load time: %d ms.\n ", (int) (nLoadTime * 1000));
    if ((Channels != 0) && (Width != 0) && (Height != 0)) {
        float *outputImg = (float *) stbi__malloc(Width * Channels * Height * sizeof(float));
        if (inputImage == NULL || outputImg == NULL) {
            fprintf(stderr, "load: %s fail!\n ", in_file);
            return -1;
        }
        for (int i = 0; i < Width * Channels * Height; ++i) {
            outputImg[i] = inputImage[i];
        }
        float Sigma = 3.0f;
        startTime = now();
        GaussianBlur(outputImg, outputImg, Width, Height, Channels, Width * Channels, Sigma);
        double nProcessTime = calcElapsed(startTime, now());
        printf("process time: %d ms.\n ", (int) (nProcessTime * 1000));
        for (int i = 0; i < Width * Channels * Height; ++i) {
            inputImage[i] = (unsigned char) (outputImg[i]);
        }
        startTime = now();
        saveImage(out_file, Width, Height, Channels, inputImage);
        double nSaveTime = calcElapsed(startTime, now());
        printf("save time: %d ms.\n ", (int) (nSaveTime * 1000));
        stbi_image_free(outputImg);
        stbi_image_free(inputImage);
    } else {
        fprintf(stderr, "load: %s fail!\n", in_file);
    }
    printf("press any key to exit. \n");
    getchar();
    return (EXIT_SUCCESS);
}
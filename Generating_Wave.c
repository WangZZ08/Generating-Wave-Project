#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>

#define PI 3.14159265359
#define MAX_AMPLITUDE(bitDepth) ((bitDepth == 8) ? 127 : pow(2.0, (bitDepth) - 1) - 1)

// WAV 標頭寫入函數
void write_wav_header(FILE *file, int sampleRate, int bitDepth, int numChannels, int numSamples) {
    int byteRate = sampleRate * numChannels * bitDepth / 8;
    int blockAlign = numChannels * bitDepth / 8;
    int subchunk2Size = numSamples * numChannels * bitDepth / 8;
    int chunkSize = 36 + subchunk2Size;

    fwrite("RIFF", 1, 4, file);
    fwrite(&chunkSize, 4, 1, file);
    fwrite("WAVE", 1, 4, file);

    fwrite("fmt ", 1, 4, file);
    int subchunk1Size = 16;
    short audioFormat = 1; // PCM 格式
    fwrite(&subchunk1Size, 4, 1, file);
    fwrite(&audioFormat, 2, 1, file);
    fwrite(&numChannels, 2, 1, file);
    fwrite(&sampleRate, 4, 1, file);
    fwrite(&byteRate, 4, 1, file);
    fwrite(&blockAlign, 2, 1, file);
    fwrite(&bitDepth, 2, 1, file);

    fwrite("data", 1, 4, file);
    fwrite(&subchunk2Size, 4, 1, file);
}

// 均勻線性量化函數
double uniform_linear_quantize(double x, double maxAmplitude) {
    return round(x * maxAmplitude) / maxAmplitude;
}

// 各種波形的生成函數原型
double sine_wave(double t, double A, double f);
double sawtooth_wave(double t, double A, double f);
double square_wave(double t, double A, double f);
double triangle_wave(double t, double A, double f);

// 主生成波形函數
void generate_wave(double *signal, double *quantized, FILE *file, const char *wave_type, int fs, double A, double f, int T, int bitDepth, int numChannels);
double calculate_sqnr(double *signal, double *quantized, int size);

int main(int argc, char *argv[]) {
    if (argc < 8) {
        fprintf(stderr, "Usage: %s <fs> <bitDepth> <numChannels> <wave_type> <frequency> <amplitude> <duration>\n", argv[0]);
        return 1;
    }

    int fs = atoi(argv[1]);
    int bitDepth = atoi(argv[2]);
    int numChannels = atoi(argv[3]);
    const char *wave_type = argv[4];
    double frequency = atof(argv[5]);
    double amplitude = atof(argv[6]);
    int duration = atoi(argv[7]);

    int totalSamples = fs * duration;
    double *signal = malloc(totalSamples * sizeof(double));
    double *quantized = malloc(totalSamples * sizeof(double));

    FILE *outputFile = stdout;

    write_wav_header(outputFile, fs, bitDepth, numChannels, totalSamples);
    generate_wave(signal, quantized, outputFile, wave_type, fs, amplitude, frequency, duration, bitDepth, numChannels);

    double sqnr = calculate_sqnr(signal, quantized, totalSamples);
    fprintf(stderr, "SQNR: %.15f\n", sqnr);

    fclose(outputFile);
    free(signal);
    free(quantized);

    return 0;
}

/* ----------------------------------- 生成 wav -----------------------------------*/
void generate_wave(double *signal, double *quantized, FILE *file, const char *wave_type, int fs, double A, double f, int T, int bitDepth, int numChannels) {
    int totalSamples = fs * T;
    double Ts = 1.0 / fs;
    double maxAmplitude = MAX_AMPLITUDE(bitDepth);

    for (int n = 0; n < totalSamples; n++) {
        double t = n * Ts;

        if (strcmp(wave_type, "sine") == 0) {
            signal[n] = sine_wave(t, A, f);
        } 
        else if (strcmp(wave_type, "sawtooth") == 0) {
            signal[n] = sawtooth_wave(t, A, f);
        } 
        else if (strcmp(wave_type, "square") == 0) {
            signal[n] = square_wave(t, A, f);
        } 
        else if (strcmp(wave_type, "triangle") == 0) {
            signal[n] = triangle_wave(t, A, f);
        }

        // 均勻線性量化
        quantized[n] = uniform_linear_quantize(signal[n], maxAmplitude);

        // 寫入量化後的數據
        if (bitDepth == 8) {
            uint8_t sample = (uint8_t)((quantized[n] + 1.0) * 127.5);
            fwrite(&sample, sizeof(uint8_t), 1, file);
        } 
        else if (bitDepth == 16) {
            int16_t sample = (int16_t)(quantized[n] * maxAmplitude);
            fwrite(&sample, sizeof(int16_t), 1, file);
        } 
        else if (bitDepth == 32) {
            int32_t sample = (int32_t)(quantized[n] * maxAmplitude);
            fwrite(&sample, sizeof(int32_t), 1, file);
        }
    }
}

// 各種波形的生成函數實現
double sine_wave(double t, double A, double f) {
    return A * sin(2.0 * PI * f * t);
}

double sawtooth_wave(double t, double A, double f) {
    double period = 1.0 / f;
    double phase = fmod(t, period);
    double normalizedPhase = phase / period;
    double sawtoothValue = 2.0 * normalizedPhase - 1.0;
    return A * sawtoothValue;
}

double square_wave(double t, double A, double f) {
    double value = sin(2.0 * PI * f * t);
    return (value >= 0) ? A * 1.0 : A * -1.0;
}

double triangle_wave(double t, double A, double f) {
    double phase = fmod(t * f, 1.0);
    double value = (phase < 0.5) ? 4.0 * phase - 1.0 : 3.0 - 4.0 * phase;
    return A * value;
}

/* ----------------------------------- 計算 sqnr -----------------------------------*/
double calculate_sqnr(double *signal, double *quantized, int size) {
    double signalPower = 0.0;
    double noisePower = 0.0;

    for (int i = 0; i < size; i++) {
        signalPower += pow(signal[i], 2);
        noisePower += pow(signal[i] - quantized[i], 2);
    }

    return 10 * log10((signalPower / size) / (noisePower / size));
}

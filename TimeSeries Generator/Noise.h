#pragma once
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <random>
#include <functional>
#include <algorithm>
#include <numeric>

enum class DataType { FLOAT, INT };


class NOISE {
public:
    NOISE() = default;

    // Вспомогательная функция для генерации шума с заданной спектральной плотностью мощности
    template<typename Func>
    static std::vector<double> generate_noise_psd(
        int N,
        double noise_std,
        Func psd_func,
        DataType dtype = DataType::FLOAT);

    // Белый шум - СПМ ~ 1
    static std::vector<double> white_noise(int N, double noise_std, DataType dtype = DataType::FLOAT);

    // Голубой шум - СПМ ~ f
    static std::vector<double> blue_noise(int N, double noise_std, DataType dtype = DataType::FLOAT);

    // Фиолетовый шум - СПМ ~ f^2
    static std::vector<double> violet_noise(int N, double noise_std, DataType dtype = DataType::FLOAT);

    // Красный шум - СПМ ~ 1/f^2
    static std::vector<double> brownian_noise(int N, double noise_std, DataType dtype = DataType::FLOAT);

    // Розовый шум - СПМ ~ 1/f
    static std::vector<double> pink_noise(int N, double noise_std, DataType dtype = DataType::FLOAT);

    ~NOISE() = default;

};
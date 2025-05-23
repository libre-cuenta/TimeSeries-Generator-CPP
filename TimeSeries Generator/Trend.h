#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <algorithm>
#include <numeric>
#include <stdexcept>

// Примечание: для полной функциональности требуются библиотеки:
// - Eigen для работы с векторами и матрицами
// - matplotlib-cpp для построения графиков
// - alglib для кубических сплайнов
// В этом коде я реализую основную логику без этих зависимостей

enum class DataType { FLOAT, INT };

using vec = std::vector<double>;

struct SplineSet {
    double a;
    double b;
    double c;
    double d;
    double x;
};

std::vector<SplineSet> spline(const vec& x, const vec& y);

vec CubicSpline(std::vector<SplineSet>& coefficient, vec& x_values);

class TREND {
public:
    TREND() = default;

    // Полиномиальный тренд - a[0] + a[1]*x + a[2]*x^2 ...
    static std::vector<double> polinom_trend(const std::vector<double>& time,
        const std::vector<double>& a = { 1.0 },
        DataType dtype = DataType::FLOAT);

    // Экспоненциальный тренд - a * e^(b*x)
    static std::vector<double> exp_trend(const std::vector<double>& time,
        double a = 1.0,
        double b = 1.0,
        DataType dtype = DataType::FLOAT);

    // Экспоненциальный по основанию тренд - a * core^(b*x)
    static std::vector<double> exp2_trend(const std::vector<double>& time,
        double core,
        double a = 1.0,
        double b = 1.0,
        DataType dtype = DataType::FLOAT);

    // Логарифмический тренд - a + b*log(c*x)
    static std::vector<double> log_trend(const std::vector<double>& time,
        double a = 1.0,
        double b = 1.0,
        double c = 1.0,
        DataType dtype = DataType::FLOAT);

    // Степенной тренд - a * (x^b)
    static std::vector<double> extend_trend(const std::vector<double>& time,
        double a = 1.0,
        double b = 1.0,
        DataType dtype = DataType::FLOAT);

    // Сплайн по точкам
    // Примечание: для полной реализации требуется библиотека для кубических сплайнов
    static std::pair<std::vector<double>, std::vector<double>> spline_trend(
        const std::vector<double>& x_points,
        const std::vector<double>& y_points,
        int num_points = 100);

    ~TREND() = default;
};
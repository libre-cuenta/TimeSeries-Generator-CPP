#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <algorithm>
#include <complex>
#include <stdexcept>

enum class DataType { FLOAT, INT };


class SEASON {
public:
    SEASON() = default;

    // Тригонометрический ряд вида - a*|sin(b*x+c)|^d * sign(sin(b*x+c))
    static std::vector<double> trigonometric_row_1(
        const std::vector<double>& time,
        double a = 1.0,
        double b = 1.0,
        double c = 0.0,
        double d = 1.0,
        DataType dtype = DataType::FLOAT);

    // Тригонометрический ряд вида - a0 + Sum(i от 1 до N): [a[i]*cos(alpha*i*(x^delta)) + b[i]*sin(alpha*i*(x^delta))]
    static std::vector<double> trigonometric_row_2(
        const std::vector<double>& time,
        double a0 = 0.0,
        const std::vector<double>& a = { 1.0 },
        const std::vector<double>& b = { 1.0 },
        double alpha = 1.0,
        double delta = 1.0,
        DataType dtype = DataType::FLOAT);

    // Поличастотная функция sin - Sum(i от 1 до N): a[i] * sin(alpha[i]*x)
    static std::vector<double> frequency_function_sin(
        const std::vector<double>& time,
        double a0 = 0.0,
        const std::vector<double>& a = { 1.0 },
        const std::vector<double>& alpha = { 1.0 },
        DataType dtype = DataType::FLOAT);

    // Поличастотная функция cos - Sum(i от 1 до N): a[i] * cos(alpha[i]*x)
    static std::vector<double> frequency_function_cos(
        const std::vector<double>& time,
        double a0 = 0.0,
        const std::vector<double>& a = { 1.0 },
        const std::vector<double>& alpha = { 1.0 },
        DataType dtype = DataType::FLOAT);

    // Псевдопериодическая функция с изменяющейся амплитудой - f(x)*sin(b*x+c)
    static std::vector<double> variable_amplitude(
        const std::vector<double>& time,
        const std::function<double(double)>& f,
        double b = 1.0,
        double c = 0.0,
        DataType dtype = DataType::FLOAT);

    // Обобщенный ряд Фурье - Sum(i от 1 до N): c[i]*exp^(j*x*lamb[i])
    static std::vector<double> furier_row(
        const std::vector<double>& time,
        const std::vector<std::complex<double>>& c,
        const std::vector<double>& lamb,
        DataType dtype = DataType::FLOAT);

    // Модулированный сигнал вида - (a0 + sin(f*x))*sin(x)
    static std::vector<double> moduling_signal(
        const std::vector<double>& time,
        double a0 = 1.0,
        double f = 1.0,
        DataType dtype = DataType::FLOAT);

    // Модулированный сигнал вида - sin(alpha*x)*cos(beta*x)
    static std::vector<double> moduling_signal2(
        const std::vector<double>& time,
        double alpha = 1.0,
        double beta = 1.0,
        DataType dtype = DataType::FLOAT);

    // Функция Вейерштрасса - Sum(i от 1 до N): (alpha^i) * cos( (beta^i)*Pi*x )
    static std::vector<double> weierstrass(
        const std::vector<double>& time,
        int N,
        double alpha = 1.0,
        double beta = 1.0,
        DataType dtype = DataType::FLOAT);

    // Линейная частотная модуляция - a0 * cos( phi0 + 2*Pi*(f0*t + (b/2)*t^2) )
    static std::vector<double> LFM(
        const std::vector<double>& time,
        double a0 = 1.0,
        double phi0 = 0.0,
        double f0 = 1.0,
        double b = 1.0,
        DataType dtype = DataType::FLOAT);

    // Ряд Фурье
    static std::vector<double> furier_row(const std::vector<double>& time,
        const std::vector<double>& c = { 1.0 },
        const std::vector<double>& lamb = { 1.0 });

    // Модулированный сигнал
    static std::vector<double> moduling_signal(const std::vector<double>& time,
        double a0 = 1.0,
        double f = 1.0);

    // Модулированный сигнал 2
    static std::vector<double> moduling_signal2(const std::vector<double>& time,
        double alpha = 1.0,
        double beta = 1.0);

    // Функция Вейерштрасса
    static std::vector<double> weierstrass(const std::vector<double>& time,
        int N,
        double alpha = 1.0,
        double beta = 1.0);

    // Линейно-частотно модулированный сигнал (ЛЧМ)
    static std::vector<double> LFM(const std::vector<double>& time,
        double a0 = 0.0,
        double phi0 = 0.0,
        double f0 = 1.0,
        double b = 2.0);

    // Функция с косинусами и синусами
    static std::vector<double> fourier_trend(const std::vector<double>& time,
        double a0 = 1.0,
        const std::vector<double>& a = { 1.0 },
        const std::vector<double>& b = { 1.0 },
        double alpha = 1.0,
        double delta = 1.0,
        DataType dtype = DataType::FLOAT);

    // Функция с синусами
    static std::vector<double> frequency_function_sin(const std::vector<double>& time,
        double a0 = 1.0,
        const std::vector<double>& a = { 1.0 },
        const std::vector<double>& alpha = { 1.0 },
        DataType dtype = DataType::FLOAT);

    // Функция с косинусами и синусами
    static std::vector<double> fourier_trend(const std::vector<double>& time,
        double a0 = 1.0,
        const std::vector<double>& a = { 1.0 },
        const std::vector<double>& b = { 1.0 },
        double alpha = 1.0,
        double delta = 1.0,
        DataType dtype = DataType::FLOAT);

    // Функция с синусами
    static std::vector<double> frequency_function_sin(const std::vector<double>& time,
        double a0 = 1.0,
        const std::vector<double>& a = { 1.0 },
        const std::vector<double>& alpha = { 1.0 },
        DataType dtype = DataType::FLOAT);

    // Функция с переменной амплитудой
    static std::vector<double> variable_amplitude(const std::vector<double>& time,
        const std::function<double(double)>& f = [](double t) { return t; },
        double b = 1.0,
        double c = 1.0,
        DataType dtype = DataType::FLOAT);

    // Функция с косинусами
    static std::vector<double> frequency_function_cos(const std::vector<double>& time,
        double a0 = 1.0,
        const std::vector<double>& a = { 1.0 },
        const std::vector<double>& alpha = { 1.0 },
        DataType dtype = DataType::FLOAT);

    ~SEASON() = default;
};

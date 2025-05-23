#include "Trend.h"

std::vector<SplineSet> spline(const vec& x, const vec& y)
{
    int n = x.size() - 1;
    vec a;
    a.insert(a.begin(), y.begin(), y.end());
    vec b(n);
    vec d(n);
    vec h;

    for (int i = 0; i < n; ++i)
        h.push_back(x[i + 1] - x[i]);

    vec alpha;
    alpha.push_back(0);
    for (int i = 1; i < n; ++i)
        alpha.push_back(3 * (a[i + 1] - a[i]) / h[i] - 3 * (a[i] - a[i - 1]) / h[i - 1]);

    vec c(n + 1);
    vec l(n + 1);
    vec mu(n + 1);
    vec z(n + 1);
    l[0] = 1;
    mu[0] = 0;
    z[0] = 0;

    for (int i = 1; i < n; ++i)
    {
        l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    l[n] = 1;
    z[n] = 0;
    c[n] = 0;

    for (int j = n - 1; j >= 0; --j)
    {
        c[j] = z[j] - mu[j] * c[j + 1];
        b[j] = (a[j + 1] - a[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3;
        d[j] = (c[j + 1] - c[j]) / 3 / h[j];
    }

    std::vector<SplineSet> output_set(n);
    for (int i = 0; i < n; ++i)
    {
        output_set[i].a = a[i];
        output_set[i].b = b[i];
        output_set[i].c = c[i];
        output_set[i].d = d[i];
        output_set[i].x = x[i];
    }
    return output_set;
}

vec CubicSpline(std::vector<SplineSet>& coefficient, vec& x_values) {
    if (coefficient.size() != x_values.size()) {
        throw std::invalid_argument("Количество коэффициентов и значений x должно совпадать");
    }

    vec result;
    for (int i = 0; i < coefficient.size(); ++i) {
        result.push_back(coefficient[i].a + coefficient[i].b * (x_values[i] - coefficient[i].x) + coefficient[i].c * std::pow(((x_values[i] - coefficient[i].x)), 2) + coefficient[i].d * (std::pow(((x_values[i] - coefficient[i].x)), 3)));
    }

    return result;
}


// Полиномиальный тренд - a[0] + a[1]*x + a[2]*x^2 ...
std::vector<double> TREND::polinom_trend(const std::vector<double>& time,
    const std::vector<double>& a = { 1.0 },
    DataType dtype = DataType::FLOAT) {
    if (time.empty()) {
        throw std::runtime_error("Необходимо задать временной промежуток для генерации");
    }

    double value = 0.0;
    std::vector<double> Y;
    Y.reserve(time.size());

    if (a.size() == 1) {
        Y.assign(time.size(), a[0]);
    }
    else {
        for (double x : time) {
            value = 0.0;
            for (size_t i = 0; i < a.size(); ++i) {
                value += a[i] * std::pow(x, i);
            }
            Y.push_back(value);
        }
    }

    if (dtype == DataType::INT) {
        for (auto& y : Y) {
            y = static_cast<int>(y);
        }
    }

    return Y;
}

// Экспоненциальный тренд - a * e^(b*x)
std::vector<double> TREND::exp_trend(const std::vector<double>& time,
    double a = 1.0,
    double b = 1.0,
    DataType dtype = DataType::FLOAT) {
    if (time.empty()) {
        throw std::runtime_error("Необходимо задать временной промежуток для генерации");
    }

    std::vector<double> Y;
    Y.reserve(time.size());

    for (double x : time) {
        Y.push_back(a * std::exp(b * x));
    }

    if (dtype == DataType::INT) {
        for (auto& y : Y) {
            y = static_cast<int>(y);
        }
    }

    return Y;
}

// Экспоненциальный по основанию тренд - a * core^(b*x)
std::vector<double> TREND::exp2_trend(const std::vector<double>& time,
    double core,
    double a = 1.0,
    double b = 1.0,
    DataType dtype = DataType::FLOAT) {
    if (time.empty()) {
        throw std::runtime_error("Необходимо задать временной промежуток для генерации");
    }

    std::vector<double> Y;
    Y.reserve(time.size());

    for (double x : time) {
        Y.push_back(a * std::pow(core, b * x));
    }

    if (dtype == DataType::INT) {
        for (auto& y : Y) {
            y = static_cast<int>(y);
        }
    }

    return Y;
}

// Логарифмический тренд - a + b*log(c*x)
std::vector<double> TREND::log_trend(const std::vector<double>& time,
    double a = 1.0,
    double b = 1.0,
    double c = 1.0,
    DataType dtype = DataType::FLOAT) {
    if (time.empty()) {
        throw std::runtime_error("Необходимо задать временной промежуток для генерации");
    }

    std::vector<double> Y;
    Y.reserve(time.size());

    for (double x : time) {
        if (x > 0) {
            Y.push_back(a + b * std::log(c * x));
        }
    }

    if (dtype == DataType::INT) {
        for (auto& y : Y) {
            y = static_cast<int>(y);
        }
    }

    return Y;
}
// Степенной тренд - a * (x^b)
std::vector<double> TREND::extend_trend(const std::vector<double>& time,
    double a = 1.0,
    double b = 1.0,
    DataType dtype = DataType::FLOAT) {
    if (time.empty()) {
        throw std::runtime_error("Необходимо задать временной промежуток для генерации");
    }

    std::vector<double> Y;
    Y.reserve(time.size());

    for (double x : time) {
        Y.push_back(a * std::pow(x, b));
    }

    if (dtype == DataType::INT) {
        for (auto& y : Y) {
            y = static_cast<int>(y);
        }
    }

    return Y;
}

// Сплайн по точкам
// Примечание: для полной реализации требуется библиотека для кубических сплайнов
std::pair<std::vector<double>, std::vector<double>> TREND::spline_trend(
    const std::vector<double>& x_points,
    const std::vector<double>& y_points,
    int num_points = 100) {

    if (x_points.size() != y_points.size()) {
        throw std::invalid_argument("Количество x и y координат должно совпадать");
    }

    if (x_points.size() < 3) {
        throw std::invalid_argument("Для построения кубического сплайна необходимо минимум 3 точки");
    }

    std::vector<SplineSet> cs = spline(x_points, y_points);


    // Заглушка для демонстрации интерфейса
    double min_x = *std::min_element(x_points.begin(), x_points.end());
    double max_x = *std::max_element(x_points.begin(), x_points.end());

    vec x_new(num_points);
    vec y_new(num_points);

    double step = (max_x - min_x) / (num_points - 1);
    for (int i = 0; i < num_points; ++i) {
        x_new[i] = min_x + i * step;
    }

    y_new = CubicSpline(cs, x_new);
      
    return { x_new, y_new };
}


#include "Season.h"


// Тригонометрический ряд вида - a*|sin(b*x+c)|^d * sign(sin(b*x+c))
std::vector<double> SEASON::trigonometric_row_1(
    const std::vector<double>& time,
    double a = 1.0,
    double b = 1.0,
    double c = 0.0,
    double d = 1.0,
    DataType dtype = DataType::FLOAT) {

    if (time.empty()) {
        throw std::runtime_error("Необходимо задать временной промежуток для генерации");
    }

    std::vector<double> Y;
    Y.reserve(time.size());

    double sin_val;
    double sign_val;
    for (double x : time) {
        sin_val = std::sin(b * x + c);
        sign_val = (sin_val > 0) ? 1.0 : ((sin_val < 0) ? -1.0 : 0.0);
        Y.push_back(a * std::pow(std::abs(sin_val), d) * sign_val);
    }

    if (dtype == DataType::INT) {
        for (auto& y : Y) {
            y = static_cast<int>(y);
        }
    }

    return Y;
}

// Тригонометрический ряд вида - a0 + Sum(i от 1 до N): [a[i]*cos(alpha*i*(x^delta)) + b[i]*sin(alpha*i*(x^delta))]
std::vector<double> SEASON::trigonometric_row_2(
    const std::vector<double>& time,
    double a0 = 0.0,
    const std::vector<double>& a = { 1.0 },
    const std::vector<double>& b = { 1.0 },
    double alpha = 1.0,
    double delta = 1.0,
    DataType dtype = DataType::FLOAT) {

    if (time.empty()) {
        throw std::runtime_error("Необходимо задать временной промежуток для генерации");
    }

    std::vector<double> Y;
    Y.reserve(time.size());

    double x_pow_delta;
    double season;
    double x_pow_delta;
    // Случай когда a и b - скаляры
    if (a.size() == 1 && b.size() == 1) {
        for (double x : time) {
            x_pow_delta = std::pow(x, delta);
            Y.push_back(a0 + a[0] * std::cos(alpha * x_pow_delta) + b[0] * std::sin(alpha * x_pow_delta));
        }
    }
    // Общий случай
    else {
        for (double x : time) {
            season = 0.0;
            x_pow_delta = std::pow(x, delta);

            for (size_t a_i = 0; a_i < a.size(); ++a_i) {
                for (size_t b_i = 0; b_i < b.size(); ++b_i) {
                    size_t max_idx = std::max(a_i, b_i) + 1;
                    season += a[a_i] * std::cos(alpha * max_idx * x_pow_delta) +
                        b[b_i] * std::sin(alpha * max_idx * x_pow_delta);
                }
            }

            Y.push_back(a0 + season);
        }
    }

    if (dtype == DataType::INT) {
        for (auto& y : Y) {
            y = static_cast<int>(y);
        }
    }

    return Y;
}

// Поличастотная функция sin - Sum(i от 1 до N): a[i] * sin(alpha[i]*x)
std::vector<double> SEASON::frequency_function_sin(
    const std::vector<double>& time,
    double a0 = 0.0,
    const std::vector<double>& a = { 1.0 },
    const std::vector<double>& alpha = { 1.0 },
    DataType dtype = DataType::FLOAT) {

    if (time.empty()) {
        throw std::runtime_error("Необходимо задать временной промежуток для генерации");
    }

    std::vector<double> Y;
    Y.reserve(time.size());

    double value;
    for (double x : time) {
        value = a0;
        for (size_t i = 0; i < std::min(a.size(), alpha.size()); ++i) {
            value += a[i] * std::sin(alpha[i] * x);
        }
        Y.push_back(value);
    }

    if (dtype == DataType::INT) {
        for (auto& y : Y) {
            y = static_cast<int>(y);
        }
    }

    return Y;
}

// Поличастотная функция cos - Sum(i от 1 до N): a[i] * cos(alpha[i]*x)
std::vector<double> SEASON::frequency_function_cos(
    const std::vector<double>& time,
    double a0 = 0.0,
    const std::vector<double>& a = { 1.0 },
    const std::vector<double>& alpha = { 1.0 },
    DataType dtype = DataType::FLOAT) {

    if (time.empty()) {
        throw std::runtime_error("Необходимо задать временной промежуток для генерации");
    }

    std::vector<double> Y;
    Y.reserve(time.size());

    double value;
    for (double x : time) {
        value = a0;
        for (size_t i = 0; i < std::min(a.size(), alpha.size()); ++i) {
            value += a[i] * std::cos(alpha[i] * x);
        }
        Y.push_back(value);
    }

    if (dtype == DataType::INT) {
        for (auto& y : Y) {
            y = static_cast<int>(y);
        }
    }

    return Y;
}

// Псевдопериодическая функция с изменяющейся амплитудой - f(x)*sin(b*x+c)
std::vector<double> SEASON::variable_amplitude(
    const std::vector<double>& time,
    const std::function<double(double)>& f,
    double b = 1.0,
    double c = 0.0,
    DataType dtype = DataType::FLOAT) {

    if (time.empty()) {
        throw std::runtime_error("Необходимо задать временной промежуток для генерации");
    }

    std::vector<double> Y;
    Y.reserve(time.size());

    for (double x : time) {
        Y.push_back(f(x) * std::sin(b * x + c));
    }

    if (dtype == DataType::INT) {
        for (auto& y : Y) {
            y = static_cast<int>(y);
        }
    }

    return Y;
}

// Обобщенный ряд Фурье - Sum(i от 1 до N): c[i]*exp^(j*x*lamb[i])
std::vector<double> SEASON::furier_row(
    const std::vector<double>& time,
    const std::vector<std::complex<double>>& c,
    const std::vector<double>& lamb,
    DataType dtype = DataType::FLOAT) {

    if (time.empty()) {
        throw std::runtime_error("Необходимо задать временной промежуток для генерации");
    }

    std::vector<double> Y;
    Y.reserve(time.size());

    std::complex<double> sum;
    std::complex<double> j;
    for (double x : time) {
        sum = std::complex<double>(0.0, 0.0);
        for (size_t i = 0; i < std::min(c.size(), lamb.size()); ++i) {
            j = std::complex<double>(0.0, 1.0);
            sum += c[i] * std::exp(j * x * lamb[i]);
        }
        Y.push_back(sum.real());
    }

    if (dtype == DataType::INT) {
        for (auto& y : Y) {
            y = static_cast<int>(y);
        }
    }

    return Y;
}

// Модулированный сигнал вида - (a0 + sin(f*x))*sin(x)
std::vector<double> SEASON::moduling_signal(
    const std::vector<double>& time,
    double a0 = 1.0,
    double f = 1.0,
    DataType dtype = DataType::FLOAT) {

    if (time.empty()) {
        throw std::runtime_error("Необходимо задать временной промежуток для генерации");
    }

    std::vector<double> Y;
    Y.reserve(time.size());

    for (double x : time) {
        Y.push_back((a0 + std::sin(f * x)) * std::sin(x));
    }

    if (dtype == DataType::INT) {
        for (auto& y : Y) {
            y = static_cast<int>(y);
        }
    }

    return Y;
}

// Модулированный сигнал вида - sin(alpha*x)*cos(beta*x)
std::vector<double> SEASON::moduling_signal2(
    const std::vector<double>& time,
    double alpha = 1.0,
    double beta = 1.0,
    DataType dtype = DataType::FLOAT) {

    if (time.empty()) {
        throw std::runtime_error("Необходимо задать временной промежуток для генерации");
    }
    std::vector<double> Y;
    Y.reserve(time.size());

    for (double x : time) {
        Y.push_back(std::sin(alpha * x) * std::cos(beta * x));
    }

    if (dtype == DataType::INT) {
        for (auto& y : Y) {
            y = static_cast<int>(y);
        }
    }

    return Y;
}

// Функция Вейерштрасса - Sum(i от 1 до N): (alpha^i) * cos( (beta^i)*Pi*x )
std::vector<double> SEASON::weierstrass(
    const std::vector<double>& time,
    int N,
    double alpha = 1.0,
    double beta = 1.0,
    DataType dtype = DataType::FLOAT) {

    if (time.empty()) {
        throw std::runtime_error("Необходимо задать временной промежуток для генерации");
    }

    std::vector<double> Y;
    Y.reserve(time.size());

    const double PI = 3.14159265358979323846;

    double sum;
    for (double x : time) {
        sum = 0.0;
        for (int i = 1; i <= N; ++i) {
            sum += std::pow(alpha, i) * std::cos(std::pow(beta, i) * PI * x);
        }
        Y.push_back(sum);
    }

    if (dtype == DataType::INT) {
        for (auto& y : Y) {
            y = static_cast<int>(y);
        }
    }

    return Y;
}

// Линейная частотная модуляция - a0 * cos( phi0 + 2*Pi*(f0*t + (b/2)*t^2) )
std::vector<double> SEASON::LFM(
    const std::vector<double>& time,
    double a0 = 1.0,
    double phi0 = 0.0,
    double f0 = 1.0,
    double b = 1.0,
    DataType dtype = DataType::FLOAT) {

    if (time.empty()) {
        throw std::runtime_error("Необходимо задать временной промежуток для генерации");
    }

    std::vector<double> Y;
    Y.reserve(time.size());

    const double PI = 3.14159265358979323846;

    for (double t : time) {
        Y.push_back(a0 * std::cos(phi0 + 2.0 * PI * (f0 * t + (b / 2.0) * t * t)));
    }

    if (dtype == DataType::INT) {
        for (auto& y : Y) {
            y = static_cast<int>(y);
        }
    }

    return Y;
}

// Ряд Фурье
std::vector<double> SEASON::furier_row(const std::vector<double>& time,
    const std::vector<double>& c = { 1.0 },
    const std::vector<double>& lamb = { 1.0 }) {
    std::vector<double> Y;
    Y.reserve(time.size());

    std::complex<double> exp_val;
    // Случай когда оба параметра скаляры
    if (c.size() == 1 && lamb.size() == 1) {
        for (double x : time) {
            exp_val = std::exp(std::complex<double>(0, 1) * x * lamb[0]);
            Y.push_back(std::abs(c[0] * exp_val));
        }
        return Y;
    }

    std::complex<double> season;
    // Случай когда один из параметров имеет длину 1
    if (c.size() == 1 || lamb.size() == 1) {
        for (double x : time) {
            season = std::complex<double>(0, 0);
            for (size_t c_i = 0; c_i < c.size(); ++c_i) {
                for (size_t lamb_i = 0; lamb_i < lamb.size(); ++lamb_i) {
                    season += c[c_i] * std::exp(std::complex<double>(0, 1) * x * lamb[lamb_i]);
                }
            }
            Y.push_back(std::abs(season));
        }
        return Y;
    }

    // Проверка входных данных
    if (c.size() != lamb.size()) {
        throw std::invalid_argument("Количество коэффициентов c и lamb должно совпадать");
    }

    std::complex<double> season;
    // Общий случай
    for (double x : time) {
        season = std::complex<double>(0, 0);
        for (size_t i = 0; i < c.size(); ++i) {
            season += c[i] * std::exp(std::complex<double>(0, 1) * x * lamb[i]);
        }
        Y.push_back(std::abs(season));
    }

    return Y;
}

// Модулированный сигнал
std::vector<double> SEASON::moduling_signal(const std::vector<double>& time,
    double a0 = 1.0,
    double f = 1.0) {
    std::vector<double> Y;
    Y.reserve(time.size());

    for (double x : time) {
        Y.push_back((a0 + std::sin(f * x)) * std::sin(x));
    }

    return Y;
}

// Модулированный сигнал 2
std::vector<double> SEASON::moduling_signal2(const std::vector<double>& time,
    double alpha = 1.0,
    double beta = 1.0) {
    std::vector<double> Y;
    Y.reserve(time.size());

    for (double x : time) {
        Y.push_back(std::sin(alpha * x) * std::cos(beta * x));
    }

    return Y;
}

// Функция Вейерштрасса
std::vector<double> SEASON::weierstrass(const std::vector<double>& time,
    int N,
    double alpha = 1.0,
    double beta = 1.0) {
    std::vector<double> Y;
    Y.reserve(time.size());

    double season;
    for (double x : time) {
        season = 0.0;
        for (int i = 1; i <= N; ++i) {
            season += std::pow(alpha, i) * std::cos(std::pow(beta, i) * M_PI * x);
        }
        Y.push_back(season);
    }

    return Y;
}

// Линейно-частотно модулированный сигнал (ЛЧМ)
std::vector<double> SEASON::LFM(const std::vector<double>& time,
    double a0 = 0.0,
    double phi0 = 0.0,
    double f0 = 1.0,
    double b = 2.0) {
    std::vector<double> Y;
    Y.reserve(time.size());

    for (double x : time) {
        Y.push_back(a0 * std::cos(phi0 + 2 * M_PI * (f0 * x + (b / 2) * std::pow(x, 2))));
    }

    return Y;
}

// Функция с косинусами и синусами
std::vector<double> SEASON::fourier_trend(const std::vector<double>& time,
    double a0 = 1.0,
    const std::vector<double>& a = { 1.0 },
    const std::vector<double>& b = { 1.0 },
    double alpha = 1.0,
    double delta = 1.0,
    DataType dtype = DataType::FLOAT) {
    std::vector<double> Y;
    Y.reserve(time.size());

    // Проверка входных данных
    if (a.size() != b.size()) {
        throw std::invalid_argument("Количество x и y координат должно совпадать");
    }

    double season;
    for (double x : time) {
        season = 0.0;
        for (size_t i = 0; i < a.size(); ++i) {
            season += a[i] * std::cos(alpha * (i + 1) * std::pow(x, delta)) +
                b[i] * std::sin(alpha * (i + 1) * std::pow(x, delta));
        }
        Y.push_back(a0 + season);
    }

    if (dtype == DataType::INT) {
        for (auto& y : Y) {
            y = static_cast<int>(y);
        }
    }

    return Y;
}

// Функция с синусами
std::vector<double> SEASON::frequency_function_sin(const std::vector<double>& time,
    double a0 = 1.0,
    const std::vector<double>& a = { 1.0 },
    const std::vector<double>& alpha = { 1.0 },
    DataType dtype = DataType::FLOAT) {
    std::vector<double> Y;
    Y.reserve(time.size());

    double season;
    // Если a и alpha - скаляры
    if (a.size() == 1 && alpha.size() == 1) {
        for (double x : time) {
            Y.push_back(a0 + a[0] * std::sin(alpha[0] * x));
        }
    }
    // Если alpha - скаляр, а a - вектор
    else if (alpha.size() == 1) {
        for (double x : time) {
            season = 0.0;
            for (size_t i = 0; i < a.size(); ++i) {
                season += a[i] * std::sin(alpha[0] * x);
            }
            Y.push_back(a0 + season);
        }
    }
    // Если a - скаляр, а alpha - вектор
    else if (a.size() == 1) {
        for (double x : time) {
            season = 0.0;
            for (size_t i = 0; i < alpha.size(); ++i) {
                season += a[0] * std::sin(alpha[i] * x);
            }
            Y.push_back(a0 + season);
        }
    }
    // Если один из векторов имеет длину 1
    else if (a.size() == 1 || alpha.size() == 1) {
        for (double x : time) {
            season = 0.0;
            for (size_t i = 0; i < a.size(); ++i) {
                for (size_t j = 0; j < alpha.size(); ++j) {
                    season += a[i] * std::sin(alpha[j] * x);
                }
            }
            Y.push_back(a0 + season);
        }
    }
    // Общий случай, когда оба вектора
    else {
        // Проверка входных данных
        if (a.size() != alpha.size()) {
            throw std::invalid_argument("Количество a и alpha должно совпадать");
        }

        for (double x : time) {
            season = 0.0;
            for (size_t i = 0; i < a.size(); ++i) {
                season += a[i] * std::sin(alpha[i] * x);
            }
            Y.push_back(a0 + season);
        }
    }

    if (dtype == DataType::INT) {
        for (auto& y : Y) {
            y = static_cast<int>(y);
        }
    }

    return Y;
}

// Функция с косинусами и синусами
std::vector<double> SEASON::fourier_trend(const std::vector<double>& time,
    double a0 = 1.0,
    const std::vector<double>& a = { 1.0 },
    const std::vector<double>& b = { 1.0 },
    double alpha = 1.0,
    double delta = 1.0,
    DataType dtype = DataType::FLOAT) {
    std::vector<double> Y;
    Y.reserve(time.size());

    // Проверка входных данных
    if (a.size() != b.size()) {
        throw std::invalid_argument("Количество x и y координат должно совпадать");
    }

    double season;
    for (double x : time) {
        season = 0.0;
        for (size_t i = 0; i < a.size(); ++i) {
            season += a[i] * std::cos(alpha * (i + 1) * std::pow(x, delta)) +
                b[i] * std::sin(alpha * (i + 1) * std::pow(x, delta));
        }
        Y.push_back(a0 + season);
    }

    if (dtype == DataType::INT) {
        for (auto& y : Y) {
            y = static_cast<int>(y);
        }
    }

    return Y;
}

// Функция с синусами
std::vector<double> SEASON::frequency_function_sin(const std::vector<double>& time,
    double a0 = 1.0,
    const std::vector<double>& a = { 1.0 },
    const std::vector<double>& alpha = { 1.0 },
    DataType dtype = DataType::FLOAT) {
    std::vector<double> Y;
    Y.reserve(time.size());

    double season;
    // Если a и alpha - скаляры
    if (a.size() == 1 && alpha.size() == 1) {
        for (double x : time) {
            Y.push_back(a0 + a[0] * std::sin(alpha[0] * x));
        }
    }
    // Если alpha - скаляр, а a - вектор
    else if (alpha.size() == 1) {
        for (double x : time) {
            season = 0.0;
            for (size_t i = 0; i < a.size(); ++i) {
                season += a[i] * std::sin(alpha[0] * x);
            }
            Y.push_back(a0 + season);
        }
    }
    // Если a - скаляр, а alpha - вектор
    else if (a.size() == 1) {
        for (double x : time) {
            season = 0.0;
            for (size_t i = 0; i < alpha.size(); ++i) {
                season += a[0] * std::sin(alpha[i] * x);
            }
            Y.push_back(a0 + season);
        }
    }
    // Если один из векторов имеет длину 1
    else if (a.size() == 1 || alpha.size() == 1) {
        for (double x : time) {
            season = 0.0;
            for (size_t i = 0; i < a.size(); ++i) {
                for (size_t j = 0; j < alpha.size(); ++j) {
                    season += a[i] * std::sin(alpha[j] * x);
                }
            }
            Y.push_back(a0 + season);
        }
    }
    // Общий случай, когда оба вектора
    else {
        // Проверка входных данных
        if (a.size() != alpha.size()) {
            throw std::invalid_argument("Количество a и alpha должно совпадать");
        }

        for (double x : time) {
            season = 0.0;
            for (size_t i = 0; i < a.size(); ++i) {
                season += a[i] * std::sin(alpha[i] * x);
            }
            Y.push_back(a0 + season);
        }
    }

    if (dtype == DataType::INT) {
        for (auto& y : Y) {
            y = static_cast<int>(y);
        }
    }

    return Y;
}

// Функция с переменной амплитудой
std::vector<double> SEASON::variable_amplitude(const std::vector<double>& time,
    const std::function<double(double)>& f = [](double t) { return t; },
    double b = 1.0,
    double c = 1.0,
    DataType dtype = DataType::FLOAT) {
    std::vector<double> Y;
    Y.reserve(time.size());

    for (double x : time) {
        Y.push_back(f(x) * std::sin(b * x + c));
    }

    if (dtype == DataType::INT) {
        for (auto& y : Y) {
            y = static_cast<int>(y);
        }
    }

    return Y;
}

// Функция с косинусами
std::vector<double> SEASON::frequency_function_cos(const std::vector<double>& time,
    double a0 = 1.0,
    const std::vector<double>& a = { 1.0 },
    const std::vector<double>& alpha = { 1.0 },
    DataType dtype = DataType::FLOAT) {
    std::vector<double> Y;
    Y.reserve(time.size());

    double season;
    // Если a и alpha - скаляры
    if (a.size() == 1 && alpha.size() == 1) {
        for (double x : time) {
            Y.push_back(a0 + a[0] * std::cos(alpha[0] * x));
        }
    }
    // Если alpha - скаляр, а a - вектор
    else if (alpha.size() == 1) {
        for (double x : time) {
            season = 0.0;
            for (size_t i = 0; i < a.size(); ++i) {
                season += a[i] * std::cos(alpha[0] * x);
            }
            Y.push_back(a0 + season);
        }
    }
    // Если a - скаляр, а alpha - вектор
    else if (a.size() == 1) {
        for (double x : time) {
            season = 0.0;
            for (size_t i = 0; i < alpha.size(); ++i) {
                season += a[0] * std::cos(alpha[i] * x);
            }
            Y.push_back(a0 + season);
        }
    }
    // Если один из векторов имеет длину 1
    else if (a.size() == 1 || alpha.size() == 1) {
        for (double x : time) {
            season = 0.0;
            for (size_t i = 0; i < a.size(); ++i) {
                for (size_t j = 0; j < alpha.size(); ++j) {
                    season += a[i] * std::cos(alpha[j] * x);
                }
            }
            Y.push_back(a0 + season);
        }
    }
    // Общий случай, когда оба вектора
    else {
        // Проверка входных данных
        if (a.size() != alpha.size()) {
            throw std::invalid_argument("Количество a и alpha должно совпадать");
        }

        for (double x : time) {
            season = 0.0;
            for (size_t i = 0; i < a.size(); ++i) {
                season += a[i] * std::cos(alpha[i] * x);
            }
            Y.push_back(a0 + season);
        }
    }

    if (dtype == DataType::INT) {
        for (auto& y : Y) {
            y = static_cast<int>(y);
        }
    }

    return Y;
}


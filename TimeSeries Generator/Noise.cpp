#include "Noise.h"

// Вспомогательная функция для генерации шума с заданной спектральной плотностью мощности
template<typename Func>
static std::vector<double> NOISE::generate_noise_psd<Func>(
    int N,
    double noise_std,
    Func psd_func,
    DataType dtype) {

    // Генерация белого шума
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> dist(0.0, noise_std);

    std::vector<double> white_noise(N);
    for (int i = 0; i < N; ++i) {
        white_noise[i] = dist(gen);
    }

    // Вычисление частот для FFT
    std::vector<double> frequencies(N / 2 + 1);
    for (int i = 0; i < frequencies.size(); ++i) {
        frequencies[i] = static_cast<double>(i) / N;
    }

    // Применение FFT к белому шуму
    std::vector<std::complex<double>> X_white = fft(white_noise);

    // Вычисление спектральной плотности мощности
    std::vector<double> S(frequencies.size());
    for (int i = 0; i < S.size(); ++i) {
        S[i] = psd_func(frequencies[i]);
    }

    // Нормализация S
    double mean_squared = 0.0;
    for (double val : S) {
        mean_squared += val * val;
    }
    mean_squared /= S.size();

    double normalization_factor = 1.0 / std::sqrt(mean_squared);
    for (double& val : S) {
        val *= normalization_factor;
    }

    // Применение спектральной плотности мощности
    std::vector<std::complex<double>> X_shaped(X_white.size());
    for (int i = 0; i < X_shaped.size(); ++i) {
        X_shaped[i] = X_white[i] * S[i];
    }

    // Обратное FFT
    std::vector<double> Y = ifft(X_shaped);

    // Преобразование типа, если требуется
    if (dtype == DataType::INT) {
        for (double& val : Y) {
            val = static_cast<int>(val);
        }
    }

    return Y;
}

// Белый шум - СПМ ~ 1
std::vector<double> NOISE::white_noise(int N, double noise_std, DataType dtype = DataType::FLOAT) {
    return generate_noise_psd(N, noise_std,
        [](double f) { return 1.0; },
        dtype);
}

// Голубой шум - СПМ ~ f
std::vector<double> NOISE::blue_noise(int N, double noise_std, DataType dtype = DataType::FLOAT) {
    return generate_noise_psd(N, noise_std,
        [](double f) { return std::sqrt(f); },
        dtype);
}

// Фиолетовый шум - СПМ ~ f^2
std::vector<double> NOISE::violet_noise(int N, double noise_std, DataType dtype = DataType::FLOAT) {
    return generate_noise_psd(N, noise_std,
        [](double f) { return f; },
        dtype);
}

// Красный шум - СПМ ~ 1/f^2
std::vector<double> NOISE::brownian_noise(int N, double noise_std, DataType dtype = DataType::FLOAT) {
    return generate_noise_psd(N, noise_std,
        [](double f) { return (f == 0) ? 0.0 : 1.0 / f; },
        dtype);
}

// Розовый шум - СПМ ~ 1/f
std::vector<double> NOISE::pink_noise(int N, double noise_std, DataType dtype = DataType::FLOAT) {
    return generate_noise_psd(N, noise_std,
        [](double f) { return (f == 0) ? 0.0 : 1.0 / std::sqrt(f); },
        dtype);
}



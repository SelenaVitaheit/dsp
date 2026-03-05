using System;

int fd = 44100;
int fs = 20000;
double T = 0.02;
int N = (int)(T * fd);
double[] chastoti = {1234, 5678, 9876, 14321, 19191, 21987};

for (int t = 0; t < chastoti.Length; t++)
{
    double f = chastoti[t];

    double[] nash = new double[N];
    for (int i = 0; i < N; i++)
    {
        double time = (double)i / fd;
        nash[i] = Math.Sin(2 * Math.PI * f * time);
    }

    short[] quant = new short[N];
    for (int i = 0; i < N; i++)
    {
        double temp = nash[i] * 32767;
        if (temp > 32767) temp = 32767;
        if (temp < -32768) temp = -32768;
        quant[i] = (short)(Math.Round(temp));
    }

    double[] sig = new double[N];
    for (int i = 0; i < N; i++)
    {
        sig[i] = (double)quant[i] / 32767;
    }

    double[] step = new double[N];
    for (int i = 0; i < N; i++)
    {
        step[i] = sig[i];
    }

    double[] Re = new double[N];
    double[] Im = new double[N];
    for (int k = 0; k < N; k++)
    {
        Re[k] = 0;
        Im[k] = 0;
        for (int i = 0; i < N; i++)
        {
            double ugol = 2 * Math.PI * k * i / N;
            Re[k] += step[i] * Math.Cos(ugol);
            Im[k] -= step[i] * Math.Sin(ugol);
        }
    }

    double[] Re_copy = new double[N];
    double[] Im_copy = new double[N];
    for (int k = 0; k < N; k++)
    {
        Re_copy[k] = Re[k];
        Im_copy[k] = Im[k];
    }

    int k_srez = (fs * N / fd);
    for (int k = k_srez; k < N - k_srez; k++)
    {
        Re[k] = 0;
        Im[k] = 0;
    }

    double[] filtr = new double[N];
    for (int i = 0; i < N; i++)
    {
        filtr[i] = 0;
        for (int k = 0; k < N; k++)
        {
            double ugol = 2 * Math.PI * k * i / N;
            filtr[i] += Re[k] * Math.Cos(ugol) - Im[k] * Math.Sin(ugol);
        }
        filtr[i] /= N;
    }

    double[] real_filter = new double[N];
    double RC = 1.0 / (2 * Math.PI * fs);
    double dt = 1.0 / fd;
    double alpha = dt / (RC + dt);
    real_filter[0] = step[0];
    for (int i = 1; i < N; i++)
    {
        real_filter[i] = real_filter[i - 1] + alpha * (step[i] - real_filter[i - 1]);
    }

    double nash_sum = 0;
    for (int i = 0; i < N; i++)
    {
        nash_sum += nash[i] * nash[i];
    }

    double ideal_err = 0;
    for (int i = 0; i < N; i++)
    {
        double err = nash[i] - filtr[i];
        ideal_err += err * err;
    }
    double ideal_isk = 100 * Math.Sqrt(ideal_err / nash_sum);

    double real_err = 0;
    for (int i = 0; i < N; i++)
    {
        double err = nash[i] - real_filter[i];
        real_err += err * err;
    }
    double real_isk = 100 * Math.Sqrt(real_err / nash_sum);
    double raznica = Math.Abs(real_isk - ideal_isk);

    double quant_sum = 0;
    for (int i = 0; i < N; i++)
    {
        double err = nash[i] - sig[i];
        quant_sum += err * err;
    }
    double quant_noise = 100 * Math.Sqrt(quant_sum / nash_sum);

    double E_before = 0;
    double E_after = 0;
    for (int k = 0; k < N; k++)
    {
        E_before += Re_copy[k] * Re_copy[k] + Im_copy[k] * Im_copy[k];
        E_after += Re[k] * Re[k] + Im[k] * Im[k];
    }

    double loss = 100 * (1 - E_after / E_before);
    if (loss < 0.001) loss = 0;

    double snr_ideal = 10 * Math.Log10(nash_sum / ideal_err);
    double snr_real = 10 * Math.Log10(nash_sum / real_err);

    Console.WriteLine("Идеальный фильтр испортил сигнал на " + ideal_isk.ToString("F3") + "%");
    Console.WriteLine("Обычный фильтр (RC) испортил на " + real_isk.ToString("F3") + "%");
    Console.WriteLine("Разница в качестве: " + raznica.ToString("F3") + "%");
    Console.WriteLine("Потеряно из-за округления: " + quant_noise.ToString("F3") + "%");
    Console.WriteLine("Отрезано высоких частот: " + loss.ToString("F3") + "%");
    Console.WriteLine("Идеальный фильтр: " + snr_ideal.ToString("F2") + " дБ (чем больше, тем лучше)");
    Console.WriteLine("Обычный фильтр: " + snr_real.ToString("F2") + " дБ");
    Console.WriteLine("");
}
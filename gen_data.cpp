#include <fstream>
#include <random>
#include <iostream>

float get_random()
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<> dis(0, 99.9); // range [0, 95)
    return dis(gen);
}

float get_random(int lower, int upper)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(lower, upper); // range [lower, upper)
    return dis(gen);
}

int main(int argc, char *argv[])
{
    if (argc < 3) 
    {
        std::cout << "./gen_data [num_rectangle] [file name]" << '\n';
        return 0;
    }
    int num_rectangle = atoi(argv[1]);

    std::fstream fout(argv[2], std::ios::out);

    float low1, low2;
    float high1, high2;

    for (int i = 0; i < num_rectangle; ++i)
    {
        low1 = get_random();
        low2 = get_random();

        high1 =  low1 + 0.01;
        high2 =  low2 + 0.01;

        // 1 means insertion, i means id
        fout << "1 " << i << " " << low1 << " " << low2 << " " << high1 << " " << high2 << '\n';
    }
    fout.close();
}
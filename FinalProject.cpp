#include <cstdlib>
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip>
#include <random>
#include <array>
#include <algorithm>
#include <functional>
#include <cmath>

using std::cout;
using std::cin;
using std::endl;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::string;

double PI = acos(-1);
string specclass;
double mass, L_sun, L_star, min_R, max_R, m_planet, mu, e, a, t_max, delt, l;

double r_max(double l_star, double l_sun)
{
    return (pow(((l_star/l_sun)/0.346), 0.5)); //AU
}

double r_min(double l_star, double l_sun)
{
    return (pow(((l_star/l_sun)/1.02), 0.5)); //AU
}

double Luminosity(double mass_star, double lumin_sun)
{
    return (pow(mass_star, 3.5) * lumin_sun); //W
}

double x(double r, double phi)
{
    return r * cos(phi); //phi in Radians
}

double y(double r, double phi)
{
    return r * sin(phi); //phi in Radians
}

int main(int argc, char** argv)
{
    m_planet = (3.0027 * pow(10, -6)); //In units of solar masses
    L_sun = 3.846 * pow(10, 26); //Watts -> Encyclopedia Britannica
    double G = 4 * pow(PI, 2);//AU^3/M_solar*yr^2
    if (argc != 2)
    {
        cout << "Please input a spectral class from the Harvard Spectral Class and rerun." << endl;
        return EXIT_SUCCESS;
    }
    specclass = argv[1];
    
    std::random_device rd;
    std::array<int, std::mt19937::state_size> seed_data;
    std::generate_n(seed_data.data(), seed_data.size(), std::ref(rd));
    std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
    std::mt19937 gen(seq);
    
    if (specclass == "o" || specclass == "O")
    {
        auto rand_r = std::bind(std::uniform_real_distribution<double>(0., 1.), gen);
        mass = pow((((rand_r() * -1.3 * 0.011684))+(pow(16, -1.3))), (1/-1.3));
    }
    
    if (specclass == "b" || specclass == "B")
    {
        auto rand_r = std::bind(std::uniform_real_distribution<double>(0., 1.), gen);
        mass = pow((((rand_r() * -1.3 * 0.2723))+(pow(2.1, -1.3))), (1/-1.3));
    }
    
    if (specclass == "a" || specclass == "A")
    {
        auto rand_r = std::bind(std::uniform_real_distribution<double>(0., 1.), gen);
        mass = pow((((rand_r() * -1.3 * 0.2035))+(pow(1.4, -1.3))), (1/-1.3));
    }
    
    if (specclass == "f" || specclass == "F")
    {
        auto rand_r = std::bind(std::uniform_real_distribution<double>(0., 1.), gen);
        mass = pow((((rand_r() * -1.3 * 0.2343))+(pow(1.04, -1.3))), (1/-1.3));

    }
    
    if (specclass == "g" || specclass == "G")
    {
        auto rand_r = std::bind(std::uniform_real_distribution<double>(0., 1.), gen);
        mass = pow((((rand_r() * -1.3 * 0.29712))+(pow(0.8, -1.3))), (1/-1.3));

    }
    
    if (specclass == "k" || specclass == "K")
    {
        auto rand_r = std::bind(std::uniform_real_distribution<double>(0., 1.), gen);
        mass = pow((((rand_r() * -1.3 * 1.144))+(pow(0.45, -1.3))), (1/-1.3));
    }
    
    if (specclass == "m" || specclass == "M")
    {
        auto rand_r = std::bind(std::uniform_real_distribution<double>(0., 1.), gen);
        mass = pow((((rand_r() * -0.3 * 2.8757))+(pow(0.08, -0.3))), (1/-0.3));
    }
    
    if (specclass != "o" && specclass != "O" && specclass != "b" && specclass != "B" && specclass != "a" && specclass != "A" && specclass != "f" && specclass != "F" && specclass != "g" && specclass != "G" && specclass != "k" && specclass != "K" && specclass != "m" && specclass != "M")
    {
        cout << "That isn't an acceptable spectral class. Please input a spectral class from the Harvard Spectral Class and rerun." << endl;
        return EXIT_SUCCESS;
    }
    
    auto rand_ecc = std::bind(std::uniform_real_distribution<double>(0, 0.3),gen);
    e = rand_ecc();
    cout << "The eccentricity is: " << e << endl;

    
    cout << "The mass of the star is: " << mass << " solar masses" << endl;
    mu = ((m_planet * mass) / (m_planet + mass));
    L_star = Luminosity(mass, L_sun);
    min_R = r_min(L_star, L_sun);
    max_R = r_max(L_star, L_sun);

    double roft_max = (min_R * ((1+e)/(1-e)));
    a = ((min_R * (1+e))/(1-pow(e, 2)));
    t_max = (sqrt(4 * pow(PI, 2) * pow(a, 3)/(G * mass)));
    delt = t_max/1000;
    cout << "The period is: " << t_max << " years" << endl;
    
    if (roft_max > max_R)
    {
        cout << "Oops! Your planet left the habitable zone :(" << endl;
        return EXIT_SUCCESS;
    }
    
    double steps = (t_max/delt);
    vector<double> roft(steps);
    vector<double> qoft(steps);
    vector<double> phioft(steps);
    
    
    roft.at(0) = min_R;
    qoft.at(0) = sqrt(mu * ((2/min_R) - (1/a)));
    phioft.at(0) = 0;

    l = sqrt((min_R * mu  * G * mass  * m_planet * (1+e)));
    //cout << l << endl;
    
    for (int i = 1; i < (steps/2); i++)
    {
        qoft.at(i) = qoft.at(i-1) - (delt * ((-(G * mass * m_planet)/(mu * pow(roft.at(i-1), 2))) + (pow(l, 2)/( pow(mu, 2) * pow(roft.at(i-1), 3)))));
        roft.at(i) = roft.at(i-1) + (delt * fabs(qoft.at(i-1)));
        phioft.at(i) = phioft.at(i-1) + (delt *  (l/(mu * pow(roft.at(i-1), 2))));
    }

    for (int i = (steps/2); i < steps; i++)
    {
        qoft.at(i) = qoft.at(i-1) + (delt * ((-(G * mass * m_planet)/(mu * pow(roft.at(i-1), 2))) + (pow(l, 2)/( pow(mu, 2) * pow(roft.at(i-1), 3)))));
        roft.at(i) = roft.at(i-1) - (delt * fabs(qoft.at(i-1)));
        phioft.at(i) = phioft.at(i-1) + (delt *  (l/(mu * pow(roft.at(i-1), 2))));
    }

    
    vector<double> xoft(roft.size());
    vector<double> yoft(roft.size());
    vector<double> t(roft.size());
    
    for (int i = 0; i<roft.size(); i++)
    {
        xoft.at(i) = roft.at(i) * cos(phioft.at(i));
        yoft.at(i) = roft.at(i) * sin(phioft.at(i));
        t.at(i) = delt * i;
    }
    
    ofstream dataout;
    string filename = "Final_Positions.txt"; //so you can change it easily
    
    //Opening a file
    dataout.open(filename);
    
    //Checking if the file can be opened
    if (!dataout.is_open())
    {
        cout << "Could not open this file; check write permission!" << endl;
        return EXIT_FAILURE;
    }
    
    //Writing to the file
    dataout << "x   y" << endl; 
    for (unsigned int i=0; i<roft.size(); i++)
    {
        dataout << xoft.at(i) << "   " << yoft.at(i) << endl;// << "   " << t.at(i) << endl;
    }
    dataout.close();
        
    return EXIT_SUCCESS;
}
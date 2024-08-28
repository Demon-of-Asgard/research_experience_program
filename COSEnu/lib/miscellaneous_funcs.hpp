// #include <math.h>
/*---------------------------------------------------------------------------*/

bool file_exists(std::string path)
{
    std::ifstream f(path.c_str(), std::ios::in);
    if (!f)
    {
        return false;
    }
    else
    {
        f.close();
        return true;
    }
}

/*---------------------------------------------------------------------------*/

float roundoff(float value, unsigned char prec)
{
  float pow_10 = pow(10.0f, (float)prec);
  return round(value * pow_10) / pow_10;
}

/*---------------------------------------------------------------------------*/

std::string draw(const int l, const string c)
{
    // write a string l times on to std::out
    std::string s = "";
    for (int i = 0; i < l; i++)
    {
        s += c;
    }
    return s;
}
/*---------------------------------------------------------------------------*/

template <typename T>
int sign(T val)
{
    return (T(0) < val) - (val < T(0));
}
/*---------------------------------------------------------------------------*/

inline void swap(FieldVar **a, FieldVar **b)
{
    FieldVar *tmp = *a;
    *a = *b;
    *b = tmp;
}

/*---------------------------------------------------------------------------*/

double g(double v, double v0, double sigma)
{
    double exponant = (v - v0) * (v - v0) / (2.0 * sigma * sigma);
    double N = sigma * sqrt(M_PI / 2.0) * (erf((1.0 + v0) / sigma / sqrt(2.0)) + erf((1.0 - v0) / sigma / sqrt(2.0)));
    return exp(-exponant) / N;
}

////////JW ADD/////////////////////////////////////////////////////////
inline double eigr(int j, int i)
{
    if(i==1){ //vp
            std::string fname = "vp_eig.txt";
        bool is_loading_0 = false;
        if (!file_exists(fname))
        {
            std::cout << "Unable to find " << fname << "\n";
            exit(1);
        }
        // std::cout << "Reading from " << fname << "\n";
        // std::printf("從文件中讀取的浮點數組 12345：\n");
        std::ifstream inputFile(fname, std::ios::in);
        if (!inputFile) {
            std::cerr << "Unable to open, JW" << std::endl;
            // return 1;
        }
        std::vector<float> numbers;
        float number;
        while (inputFile >> number) {
            numbers.push_back(number);
        }
        inputFile.close();
        return numbers[j];
    }
    if(i==0){ //vn
            std::string fname = "vn_eig.txt";
        bool is_loading_0 = false;
        if (!file_exists(fname))
        {
            std::cout << "Unable to find " << fname << "\n";
            exit(1);
        }
        // std::cout << "Reading from " << fname << "\n";
        // std::printf("從文件中讀取的浮點數組 12345：\n");
        std::ifstream inputFile(fname, std::ios::in);
        if (!inputFile) {
            std::cerr << "Unable to open, JW" << std::endl;
            // return 1;
        }
        std::vector<float> numbers;
        float number;
        while (inputFile >> number) {
            numbers.push_back(number);
        }
        inputFile.close();
        return numbers[j];
    }
    
}

////////////////////////////////////////////////////////////////////////

/*---------------------------------------------------------------------------*/

inline double eps(double z, double z0, double amp)
{
    //return amp * exp(-(z - z0) * (z - z0) / 50.0);
    // return amp * exp(-(z - z0) * (z - z0) / 0.001);

    // return amp * sin(500.0 * z);
    // return amp * (double)rand() / RAND_MAX;
    double km = 3.0;
    double Zm=20000.0;
    double rtn=0.0;
    unsigned int seed = 444;
    srand(seed);
    for(int i=0;i<20000;i++){
        rtn += cos((i+1.0) / Zm * km * z + (double)rand() / RAND_MAX * 2 * M_PI);
    }
    // return amp * (sin(0.1 * z)+sin(0.3 * z)+sin(0.6 * z)) / 3.0;
    return amp * rtn / Zm;
    // return 
    //return amp;
}
inline double epsr(double z, double z0, double amp,double rdph[20000])
{
    //return amp * exp(-(z - z0) * (z - z0) / 50.0);
    // return amp * exp(-(z - z0) * (z - z0) / 0.001);

    // return amp * sin(500.0 * z);
    // return amp * (double)rand() / RAND_MAX;
    double km = 3.0;
    double Zm=20000.0;
    double rtn=0.0;
    for(int i=0;i<(int)Zm;i++){
        // rtn += cos((i+1.0) / Zm * km * z + rdph[i]);
        rtn += cos((i+1.0) / Zm * km * z) / Zm /((i+1.0) / Zm * km);
        // rtn += cos((i+1.0) / Zm * km * z) / Zm;
        // printf("%f\n",rdph[i]);
    }
    // return amp * (sin(0.1 * z)+sin(0.3 * z)+sin(0.6 * z)) / 3.0;
    return amp * rtn;
    // return 
    //return amp;
}
inline double epsi(double z, double z0, double amp,double rdph[20000])
{
    //return amp * exp(-(z - z0) * (z - z0) / 50.0);
    // return amp * exp(-(z - z0) * (z - z0) / 0.001);

    // return amp * sin(500.0 * z);
    // return amp * (double)rand() / RAND_MAX;
    double km = 3.0;
    double Zm=20000.0;
    double rtn=0.0;
    for(int i=0;i<(int)Zm;i++){
        // rtn += sin((i+1.0) / Zm * km * z + rdph[i]);
        // rtn += sin((i+1.0) / Zm * km * z)  / Zm;
        rtn += sin((i+1.0) / Zm * km * z)  / Zm / ((i+1.0)/ Zm * km);
        // rtn += sin((i+1.0) / Zm * km * z)  / Zm;
        // printf("%f\n",rdph[i]);
    }
    // return amp * (sin(0.1 * z)+sin(0.3 * z)+sin(0.6 * z)) / 3.0;
    return amp * rtn;
    // return 
    //return amp;
}

/*---------------------------------------------------------------------------*/

inline double eps_(double z, double z0, double amp)
{
    double e = eps(z, z0, amp);
    return sqrt(1.0 - e * e);
}

/*---------------------------------------------------------------------------*/

inline double random_amp(double a)
{
    return a * (double)rand() / RAND_MAX;
}
/*---------------------------------------------------------------------------*/

double hev(double x, double x0)
{
    /* Heviside theta function */
    return 0.5 * (1.0 + sign(x - x0));
}

/*---------------------------------------------------------------------------*/

double gauss(double x, double x0, double sigma)
{
    double exponant = (x - x0) * (x - x0) / (2.0 * sigma * sigma);
    return exp(-exponant);
}

/*---------------------------------------------------------------------------*/

double L(double x, unsigned int order)
{
    // Legendre polynomials.
    double value;
    switch (order)
    {
    case 0:
        value = 1;
        break;
    case 1:
        value = x;
        break;
    case 2:
        value = 0.5 * (3 * x * x - 1);
        break;
    case 3:
        value = 0.5 * (5 * pow(x, 3) - 3 * x);
        break;
    default:
        value = 0;
    }
    return (value);
}

/*---------------------------------------------------------------------------*/


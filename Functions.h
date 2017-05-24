#include <string>

// Empirical fit to secondary electron emission equation as in Stangeby
double sec(double Te, char material);

void WarnOnce(bool &T, std::string Message);

double maxwellian(double E, double T);

// Ion backscattering
double ionback(double E, char isotope, char material, int flag);

double backscatter(double Te, double Ti, double mi, double Vion, char material, double &RE, double &RN);


double round_to_digits(double value, int digits);


/** @file Functions.h
 *  @brief File defining functions which are used globally
 *  
 *  Funciton definitions which are used across the entire scope of the project
 *  
 *  @author Luke Simons (ls5115@ic.ac.uk)
 *  @bug Functions should be placed within a namespace
 *  @bug sec function should be placed within Model
 *  @bug backscatter function should be placed within Model
 *  @bug maxwellian function should be placed within Model
 *  @bug ionback function should be placed within Model
 *  @bug ionback function should be placed within Model
 */

#ifndef __FUNCTIONS_H_INCLUDED__
#define __FUNCTIONS_H_INCLUDED__

#include <string>
#include <cmath>
#include <math.h>
#include <stdio.h>
#include <exception>

#include "Constants.h"
#include <random>

const double LOW = -PI/2;
const double HIGH = PI/2;

/** @brief Cumulative Distribution Function of cos^2(x)
 *  
 *  Function to return the cumulative distribution of cos^2(x)
 *  @param x The variable of the function cos^2(x)
 *  @return the cumulative distribution function of cos^2(x)
 */
double cos2_cdf(double x);

/** @brief Inverse Cumulative Distribution Function of cos^2(x)
 *  
 *  Function to return the inverse cumulative distribution of cos^2(x)
 *  @param u The variable of the function cos^2(x)
 *  @return the inverse cumulative distribution function of cos^2(x)
 */
double inverse_cos2_cdf(double u);

/** @brief A number randomly distributed over cos^2(x)
 *  
 *  Function to return a cos^2(x) distribution
 *  @param rng reference to the mt19937 random number generator
 *  @return a random number following from a cos^2(x) distribution
 */
double custom_distribution(std::mt19937& rng);

/** @brief Warning Message to be printed only once
 *  
 *  Function to print warning message to the screen only once. This is 
 *  achieved by passing a static bool variable which is set to true. By default
 *  the \p Message is preceeded by a warning string.
 *  @param MessageNotDisplayed is true if the message has yet to be displayed
 *  @param Message is the message to be displayed,
 */
void WarnOnce(bool &MessageNotDisplayed, std::string Message);

//!< Get the sign of a number
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

// Empirical fit to secondary electron emission equation as in Stangeby
double sec(double Te, char material);


double maxwellian(double E, double T);

// Ion backscattering
double ionback(double E, char isotope, char material, int flag);

double backscatter(double Te, double Ti, double mi, double Vion, char material, 
    double &RE, double &RN);

/** @brief Round a number to a specific number of digits
 *  
 *  Function to round \p value to a specific number of \p digits
 *  achieved by passing a static bool variable which is set to true. By default
 *  the \p Message is preceeded by a warning string.
 *  @param value the number which is to be rounded
 *  @param digits the number of digits to round to
 *  @return the number rounded to \p digits
 */
double round_to_digits(double value, int digits);

/** @brief Lambert W function. 
 *  Was ~/C/LambertW.c written K M Briggs Keith dot Briggs at bt dot com 97 
 *  May 21.  
 *  Revised KMB 97 Nov 20; 98 Feb 11, Nov 24, Dec 28; 99 Jan 13; 00 Feb 23; 
 *  01 Apr 09
 *
 *  Computes Lambert W function, principal branch.
 *  See LambertW1.c for -1 branch.
 *
 *  Returned value W(z) satisfies W(z)*exp(W(z))=z
 *  test data...
 *     W(1)= 0.5671432904097838730
 *     W(2)= 0.8526055020137254914
 *     W(20)=2.2050032780240599705
 *  To solve (a+b*R)*exp(-c*R)-d=0 for R, use
 *  R=-(b*W(-exp(-a*c/b)/b*d*c)+a*c)/b/c
 *
 *  Test: 
 *    gcc -DTESTW LambertW.c -o LambertW -lm && LambertW
 *  Library:
 *    gcc -O3 -c LambertW.c 
 *  @param z argument to LambertW funciton
 */
double LambertW(const double z);

/** @brief Lambert W Failer. 
 *  
 *  Used to throw an exception when the LambertW function fails
 *  @see LambertW()
 */
struct LambertWFailure : public std::exception{
  const char * what () const throw ();
};

#endif /* __FUNCTIONS_H_INCLUDED__ */

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include "roots.hpp"
#include <functional>
#include <algorithm>

bool bisection(std::function<double(double)> poly, double left, double right, double* root) {
    double a = std::min(left,right);
    double b = std::max(left,right);
    if ((!a)||(!b)) {
        return false;
    }
    double i = a;
    for (i = a; !(i >= b); i += 0.01) {
        if (poly(i) < 0) {
            a = i;
            break;
        }
    }
    for (i = b; !(i <= a); i -= 0.01) {
        if (poly(i) > 0) {
            b = i;
            break;
        }
    }

    // Now, a and b will have a negative and positive evaluation respectively.
    // If no negative or no positive value exists between the given boundaries,
    // The function will now realize this and return false. 

    if ((poly(a) > 0)||(poly(b) < 0)) {
        return false; 
    }

    double c = 50;
    int iterations = 0;
    while (iterations < 1000000) {
        c = (a + b)/2;
        double fc = poly(c);
        if (!((0.000001 < fc)||(fc < -0.000001))) {
            *root = c;
            return true;
        }
        if (fc > 0)
            b = c;
        else
            a = c;
        iterations += 1;
    }
    return false;
}

bool regula_falsi(std::function<double(double)> poly, double left, double right, double* root) {
    // Reuse code from bisection part.
    double a = std::min(left,right);
    double b = std::max(left,right);
    if ((!a)||(!b)) {
        return false;
    }
    double i = a;
    for (i = a; !(i >= b); i += 0.01) {
        if (poly(i) < 0) {
            a = i;
            break;
        }
    }
    for (i = b; !(i <= a); i -= 0.01) {
        if (poly(i) > 0) {
            b = i;
            break;
        }
    }

    // Now, a and b will have a negative and positive evaluation respectively.
    // If no negative or no positive value exists between the given boundaries,
    // The function will now realize this and return false. 
    if ((poly(a) > 0)||(poly(b) < 0)) {
        return false; 
    }
    double c = 0;
    int iterations = 0;
    while (iterations < 1000000) {
        c = ((a*poly(b) - b*poly(a))/(poly(b)-poly(a)));
        double fc = poly(c);
        if (!((0.000001 < fc)||(fc < -0.000001))) {
            *root = c;
            return true;
        }
        if (poly(a) * fc < 0)
            b = c;
        else
            a = c;
        iterations += 1;
    }

    return false;
}

bool newton_raphson(std::function<double(double)> poly, std::function<double(double)> dpoly, double left, double right, double guess, double* root) {
    double c = guess;
    double cnew;
    int iterations = 0;
    while (iterations < 1000000) {
        cnew = (c - poly(c)/dpoly(c));
        if ((std::abs(poly(cnew)-poly(c))<0.00000000001)&&(std::abs(dpoly(c))<0.0000000001)){
            return false;
        }
        if (std::abs(cnew - c) < 0.000001) {
            if ((left <= cnew)&&(right >= cnew)) {
                *root = cnew;
                return true;
            }
            else {
                return false;
            }
        }
        c = cnew;
        iterations += 1;
    }
    return false;
}

bool secant(std::function<double(double)> poly, double left, double right, double guess, double* root) {
    double g1 = guess;
    double g0 = (left+right)/2;
    // Test to see if its just a constant value. Technically, this is just guesswork.
    if ((poly(g1)==poly(g0))&&(poly(g1)==poly(((0.5*left)+1.5*right)/(2)))) {
        return false;
    }
    int iterations = 0;
    double g2;
    while (iterations < 1000000) {
        g2 = g1 - poly(g1)*((g1-g0)/(poly(g1)-poly(g0)));
        if (std::abs((g2 - g1)) < 0.000001) {
            if ((left <= g2)&&(right >= g2)){
                *root = g2;
                return true;
            }
            else {
                return false;
            }
        }
        if ((poly(g2)==poly(g1))&&(poly(g1)==poly(g0))) {
            return false;
        }
        g0 = g1;
        g1 = g2;
        iterations += 1;
    }
    return false;
}

int main() {
    return 0;
}
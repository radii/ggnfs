
#include <math.h>
double rint(double);
float rintf(float);

double rint(double x) {
  return ( (x<0.)? -floor(-x+.5):floor(x+.5) );
}

float rintf(float x) {
  return ( (x<(float)0.)? -(float)floor(-x+.5):(float)floor(x+.5) );
}

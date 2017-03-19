#include "dda_kinematics.h"

/** \file G-code axis system to stepper motor axis system conversion.
*/

#include <stdlib.h>

#include "dda_maths.h"
#include "sersendf.h"


void
carthesian_to_carthesian(const TARGET *startpoint, const TARGET *target,
                         axes_uint32_t delta_um, axes_int32_t steps) {
  enum axis_e i;

  for (i = X; i < E; i++) {
    delta_um[i] = (uint32_t)labs(target->axis[i] - startpoint->axis[i]);
    steps[i] = um_to_steps(target->axis[i], i);
  }

  /* Replacing the above five lines with this costs about 200 bytes binary
     size on AVR, but also takes about 120 clock cycles less during movement
     preparation. The smaller version was kept for our Arduino Nano friends.
    delta_um[X] = (uint32_t)labs(target->axis[X] - startpoint->axis[X]);
    steps[X] = um_to_steps(target->axis[X], X);
    delta_um[Y] = (uint32_t)labs(target->axis[Y] - startpoint->axis[Y]);
    steps[Y] = um_to_steps(target->axis[Y],
    Y);
    delta_um[Z] = (uint32_t)labs(target->axis[Z] - startpoint->axis[Z]);
    steps[Z] = um_to_steps(target->axis[Z], Z);
  */
  /*
    sersendf_P(PSTR("X           %lu\n"), target->axis[X]);
    sersendf_P(PSTR("Y           %lu\n"), target->axis[Y]);
    sersendf_P(PSTR("X steps     %lu\n"), steps[X]);
    sersendf_P(PSTR("Y steps     %lu\n"), steps[Y]);
    sersendf_P(PSTR("X delta_um  %lu\n"), delta_um[X]);
    sersendf_P(PSTR("Y delta_um  %lu\n"), delta_um[Y]);
  */
}
void
carthesian_to_corexy(const TARGET * startpoint, const TARGET * target,
                     axes_uint32_t delta_um, axes_int32_t steps) {

  delta_um[X] = (uint32_t)labs((target->axis[X] - startpoint->axis[X]) +
                               (target->axis[Y] - startpoint->axis[Y]));
  delta_um[Y] = (uint32_t)labs((target->axis[X] - startpoint->axis[X]) -
                               (target->axis[Y] - startpoint->axis[Y]));
  delta_um[Z] = (uint32_t)labs(target->axis[Z] - startpoint->axis[Z]);
  axes_um_to_steps_corexy(target->axis, steps);
}

void carthesian_to_cylindrical(const TARGET * startpoint, const TARGET * target,
                               axes_uint32_t delta_um, axes_int32_t steps) {
  //Let X be the radius and Y be the rotation angle theta
  //Radius sqrt(X^2 + Y^2). The approx_distance function does this already for us efficiently
  delta_um[X] = (uint32_t)labs((radius_approx(target->axis[X], target->axis[Y])) -
                               (radius_approx(startpoint->axis[X], startpoint->axis[Y])));
  //Theta = atan(Y/X) solved with CORDIC algorithm, integer computation
  //Result in millidegrees of angle
  delta_um[Y] = (uint32_t)labs((atan_CORDIC(target->axis[X], target->axis[Y])) -
                               (atan_CORDIC(startpoint->axis[X], startpoint->axis[Y])));
  delta_um[Z] = (uint32_t)labs(target->axis[Z] - startpoint->axis[Z]);

  axes_um_to_steps_cylindrical(target->axis, steps);
  /*sersendf_P(PSTR("\nX delta_um  %lu\n"), delta_um[X]);
  sersendf_P(PSTR("Y delta_md  %lu\n"), delta_um[Y]);
  sersendf_P(PSTR("X0          %lu\n"), startpoint->axis[X]);
  sersendf_P(PSTR("Y0          %lu\n"), startpoint->axis[Y]);
  sersendf_P(PSTR("X1          %lu\n"), target->axis[X]);
  sersendf_P(PSTR("Y1          %lu\n"), target->axis[Y]);*/
}

void axes_um_to_steps_cylindrical(const axes_int32_t um, axes_int32_t steps) {
  //sersendf_P(PSTR("X um        %lu\n"), um[X]);
  //sersendf_P(PSTR("Y um        %lu\n"), um[Y]);
  steps[Z] = um_to_steps(um[Z], Z);
  steps[X] = radius_approx(um[X], um[Y]);
  //sersendf_P(PSTR("radius      %lu\n"), radius_approx(um[X], um[Y]));
  steps[X] = um_to_steps(steps[X], X);
  //sersendf_P(PSTR("STEPX      %lu\n"), steps[X]);  
  steps[Y] =  atan_CORDIC(um[X], um[Y]);
  //sersendf_P(PSTR("angle       %lu\n"), steps[Y]);
  steps[Y] = um_to_steps(steps[Y], Y);
}

void axes_um_to_steps_cartesian(const axes_int32_t um, axes_int32_t steps) {
  enum axis_e i;

  for (i = X; i < E; i++) {
    steps[i] = um_to_steps(um[i], i);
  }
}

void axes_um_to_steps_corexy(const axes_int32_t um, axes_int32_t steps) {
  steps[X] = um_to_steps(um[X] + um[Y], X);
  steps[Y] = um_to_steps(um[X] - um[Y], Y);
  steps[Z] = um_to_steps(um[Z], Z);
}



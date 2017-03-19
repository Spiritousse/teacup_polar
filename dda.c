#include	"dda.h"

/** \file
	\brief Digital differential analyser - this is where we figure out which steppers need to move, and when they need to move
*/

#include	<string.h>
#include	<stdlib.h>
#include	<math.h>

#include	"dda_maths.h"
#include "preprocessor_math.h"
#include "dda_kinematics.h"
#include	"dda_lookahead.h"
#include "cpu.h"
#include	"timer.h"
#include	"serial.h"
#include	"gcode_parse.h"
#include	"dda_queue.h"
#include	"debug.h"
#include	"sersendf.h"
#include	"pinio.h"
#include "memory_barrier.h"
//#include "graycode.c"

#ifdef	DC_EXTRUDER
#include	"heater.h"
#endif


/*
	position tracking
*/

/// \var startpoint
/// \brief target position of last move in queue
TARGET BSS startpoint;

/// \var startpoint_steps
/// \brief target position of last move in queue, expressed in steps
TARGET BSS startpoint_steps;

/// \var current_position
/// \brief actual position of extruder head
/// \todo make current_position = real_position (from endstops) + offset from G28 and friends
TARGET BSS current_position;

/// \var move_state
/// \brief numbers for tracking the current state of movement
MOVE_STATE BSS move_state;

/// \var steps_per_m_P
/// \brief motor steps required to advance one meter on each axis
#ifdef KINEMATICS_CYLINDRICAL
static const axes_uint32_t PROGMEM steps_per_m_P = {
  STEPS_PER_M_X,
  STEPS_PER_ROTATION_Y,
  STEPS_PER_M_Z,
  STEPS_PER_M_E
};
#else
static const axes_uint32_t PROGMEM steps_per_m_P = {
  STEPS_PER_M_X,
  STEPS_PER_M_Y,
  STEPS_PER_M_Z,
  STEPS_PER_M_E
};
#endif

/// \var maximum_feedrate_P
/// \brief maximum allowed feedrate on each axis
static const axes_uint32_t PROGMEM maximum_feedrate_P = {
  MAXIMUM_FEEDRATE_X,
  MAXIMUM_FEEDRATE_Y,
  MAXIMUM_FEEDRATE_Z,
  MAXIMUM_FEEDRATE_E
};

/// \var c0_P
/// \brief Initialization constant for the ramping algorithm. Timer cycles for
///        first step interval.
#ifdef KINEMATICS_CYLINDRICAL
static const axes_uint32_t PROGMEM c0_P = {
  (uint32_t)((double)F_CPU / SQRT((double)STEPS_PER_M_X * ACCELERATION / 2000.)),
  (uint32_t)((double)F_CPU / SQRT((double)STEPS_PER_ROTATION_Y * ACCELERATION / 2000.)),
  (uint32_t)((double)F_CPU / SQRT((double)STEPS_PER_M_Z * ACCELERATION / 2000.)),
  (uint32_t)((double)F_CPU / SQRT((double)STEPS_PER_M_E * ACCELERATION / 2000.))
};
#else
static const axes_uint32_t PROGMEM c0_P = {
  (uint32_t)((double)F_CPU / SQRT((double)STEPS_PER_M_X * ACCELERATION / 2000.)),
  (uint32_t)((double)F_CPU / SQRT((double)STEPS_PER_M_Y * ACCELERATION / 2000.)),
  (uint32_t)((double)F_CPU / SQRT((double)STEPS_PER_M_Z * ACCELERATION / 2000.)),
  (uint32_t)((double)F_CPU / SQRT((double)STEPS_PER_M_E * ACCELERATION / 2000.))
};
#endif

/*! Set the direction of the 'n' axis
*/
static void set_direction(DDA *dda, enum axis_e n, int32_t delta) {
  uint8_t dir = (delta >= 0) ? 1 : 0;

  if (n == X)
    dda->x_direction = dir;
  else if (n == Y)
    dda->y_direction = dir;
  else if (n == Z)
    dda->z_direction = dir;
  else if (n == E)
    dda->e_direction = dir;
}

/*! Find the direction of the 'n' axis
*/
static int8_t get_direction(DDA *dda, enum axis_e n) {
  if ((n == X && dda->x_direction) ||
      (n == Y && dda->y_direction) ||
      (n == Z && dda->z_direction) ||
      (n == E && dda->e_direction))
    return 1;
  else
    return -1;
}

/*! Inititalise DDA movement structures
*/
void dda_init(void) {
  // set up default feedrate
  if (startpoint.F == 0)
    //startpoint.F = next_target.target.F = SEARCH_FEEDRATE_Z;
    startpoint.F = next_target.target.F = SEARCH_FEEDRATE_X;
  if (startpoint.e_multiplier == 0)
    startpoint.e_multiplier = next_target.target.e_multiplier = 256;
  if (startpoint.f_multiplier == 0)
    startpoint.f_multiplier = next_target.target.f_multiplier = 256;
}

/*! Distribute a new startpoint to DDA's internal structures without any movement.

	This is needed for example after homing or a G92. The new location must be in startpoint already.
*/
void dda_new_startpoint(void) {
  axes_um_to_steps(startpoint.axis, startpoint_steps.axis);
  startpoint_steps.axis[E] = um_to_steps(startpoint.axis[E], E);
  sersendf_P(PSTR("startpoint %lu\n"), startpoint.axis[X], startpoint.axis[Y]);
}


void dda_create(DDA *dda, const TARGET *target) {
  axes_uint32_t delta_um;
  axes_int32_t steps;
  uint32_t	distance, c_limit, c_limit_calc;
  enum axis_e i;
#ifdef LOOKAHEAD
  // Number the moves to identify them; allowed to overflow.
  static uint8_t idcnt = 0;
  static DDA* prev_dda = NULL;

  if ((prev_dda && prev_dda->done) || dda->waitfor_temp)
    prev_dda = NULL;
#endif

  if (dda->waitfor_temp)
    return;

  // We end at the passed target.
  memcpy(&(dda->endpoint), target, sizeof(TARGET));

  if (DEBUG_DDA && (debug_flags & DEBUG_DDA))
    sersendf_P(PSTR("\nCreate: X %lq  Y %lq  Z %lq  F %lu\n"),
               dda->endpoint.axis[X], dda->endpoint.axis[Y],
               dda->endpoint.axis[Z], dda->endpoint.F);

  // Apply feedrate multiplier.
  if (dda->endpoint.f_multiplier != 256) {
    dda->endpoint.F *= dda->endpoint.f_multiplier;
    dda->endpoint.F += 128;
    dda->endpoint.F /= 256;
  }

#ifdef LOOKAHEAD
  // Set the start and stop speeds to zero for now = full stops between
  // moves. Also fallback if lookahead calculations fail to finish in time.
  dda->crossF = 0;
  dda->start_steps = 0;
  dda->end_steps = 0;
  // Give this move an identifier.
  dda->id = idcnt++;
#endif

  // Handle bot axes. They're subject to kinematics considerations.
  code_axes_to_stepper_axes(&startpoint, target, delta_um, steps);
  //A MODIFIER POUR REMETTRE EN MARCHE AXE Z
  for (i = X; i < Z  ; i++) {
    int32_t delta_steps;

    delta_steps = steps[i] - startpoint_steps.axis[i];


    /*if (i == Y)
      {
      sersendf_P(PSTR("delta_steps Y = %ld\n"), delta_steps);
      sersendf_P(PSTR("steps[Y] = %lu\n"), steps[i]);
      sersendf_P(PSTR("steps start[Y] = %lu\n"), startpoint_steps.axis[i]);
      }*/
    dda->delta[i] = (uint32_t)labs(delta_steps);
    startpoint_steps.axis[i] = steps[i];

    set_direction(dda, i, delta_steps);

#ifdef LOOKAHEAD
    dda->delta_um[i] = (delta_steps >= 0) ?
                       (int32_t)delta_um[i] : -(int32_t)delta_um[i];
#endif
  }
  sersendf_P(PSTR("X"));
  (dda->x_direction) ? (sersendf_P(PSTR("+"))) : (sersendf_P(PSTR("-")));
  sersendf_P(PSTR("%lu"), dda->delta[X]);
  sersendf_P(PSTR("Y"));
  (dda->y_direction) ? (sersendf_P(PSTR("+"))) : (sersendf_P(PSTR("-")));
  sersendf_P(PSTR("%lu"), dda->delta[Y]);
  sersendf_P(PSTR("N\n"));

  //sersendf_P(PSTR("X dir       %u\n"), dda->x_direction);
  //sersendf_P(PSTR("Y dir       %u\n"), dda->y_direction);



  // Handle extruder axes. They act independently from the bots kinematics
  // type, but are subject to other special handling.
  steps[E] = um_to_steps(target->axis[E], E);

  // Apply extrusion multiplier.
  if (target->e_multiplier != 256) {
    steps[E] *= target->e_multiplier;
    steps[E] += 128;
    steps[E] /= 256;
  }

  if ( ! target->e_relative) {
    int32_t delta_steps;

    delta_um[E] = (uint32_t)labs(target->axis[E] - startpoint.axis[E]);
    delta_steps = steps[E] - startpoint_steps.axis[E];
    dda->delta[E] = (uint32_t)labs(delta_steps);
    startpoint_steps.axis[E] = steps[E];

    set_direction(dda, E, delta_steps);
#ifdef LOOKAHEAD
    dda->delta_um[E] = (delta_steps >= 0) ?
                       (int32_t)delta_um[E] : -(int32_t)delta_um[E];
#endif
  }
  else {
    // When we get more extruder axes:
    // for (i = E; i < AXIS_COUNT; i++) { ...
    delta_um[E] = (uint32_t)labs(target->axis[E]);
    dda->delta[E] = (uint32_t)labs(steps[E]);
#ifdef LOOKAHEAD
    dda->delta_um[E] = target->axis[E];
#endif
    dda->e_direction = (target->axis[E] >= 0) ? 1 : 0;
  }

  if (DEBUG_DDA && (debug_flags & DEBUG_DDA))
    sersendf_P(PSTR("[%ld,%ld,%ld,%ld]"),
               target->axis[X] - startpoint.axis[X], target->axis[Y] - startpoint.axis[Y],
               target->axis[Z] - startpoint.axis[Z], target->axis[E] - startpoint.axis[E]);

  for (i = X; i <  Z; i++) {
    if (i == X || dda->delta[i] > dda->total_steps) {
      dda->fast_axis = i;
      dda->total_steps = dda->delta[i];
      dda->fast_um = delta_um[i];
      dda->fast_spm = pgm_read_dword(&steps_per_m_P[i]);
      //sersendf_P(PSTR("steps %lu\n"), dda->delta[i]);
    }

  }
  //sersendf_P(PSTR("total steps %lu\n"), dda->total_steps);

  if (DEBUG_DDA && (debug_flags & DEBUG_DDA))
    sersendf_P(PSTR(" [ts:%lu"), dda->total_steps);

  if (dda->total_steps == 0) {
    dda->nullmove = 1;
  }
  else {
    // get steppers ready to go
    power_on();
    stepper_enable();
    x_enable();
    y_enable();
    z_enable();
    // #else Z is enabled in dda_start().
    e_enable();

    // since it's unusual to combine X, Y and Z changes in a single move on reprap, check if we can use simpler approximations before trying the full 3d approximation.
    if (delta_um[Z] == 0)
      distance = approx_distance(delta_um[X], delta_um[Y]);
    else if (delta_um[X] == 0 && delta_um[Y] == 0)
      distance = delta_um[Z];
    else
      distance = approx_distance_3(delta_um[X], delta_um[Y], delta_um[Z]);

    if (distance < 2)
      distance = delta_um[E];

    if (DEBUG_DDA && (debug_flags & DEBUG_DDA))
      sersendf_P(PSTR(",ds:%lu"), distance);


    // pre-calculate move speed in millimeter microseconds per step minute for less math in interrupt context
    uint32_t move_duration = ((distance * 2400) / dda->total_steps) * (F_CPU / 40000);
    //sersendf_P(PSTR("movement duration %lu \n"), move_duration);

    // similarly, find out how fast we can run our axes.
    // do this for each axis individually, as the combined speed of two or more axes can be higher than the capabilities of a single one.
    c_limit = 0;
    for (i = X; i < AXIS_COUNT; i++) {
      c_limit_calc = (delta_um[i] * 2400L) /
                     dda->total_steps * (F_CPU / 40000) /
                     pgm_read_dword(&maximum_feedrate_P[i]);
      if (c_limit_calc > c_limit)
        c_limit = c_limit_calc;
    }

    dda->c_min = move_duration / dda->endpoint.F;
    if (dda->c_min < c_limit) {
      dda->c_min = c_limit;
      dda->endpoint.F = move_duration / dda->c_min;
    }

    // Lookahead can deal with 16 bits ( = 1092 mm/s), only.
    if (dda->endpoint.F > 65535)
      dda->endpoint.F = 65535;

    // Acceleration ramps are based on the fast axis, not the combined speed.
    dda->rampup_steps =
      acc_ramp_len(muldiv(dda->fast_um, dda->endpoint.F, distance),
                   dda->fast_spm);

    if (dda->rampup_steps > dda->total_steps / 2)
      dda->rampup_steps = dda->total_steps / 2;
    dda->rampdown_steps = dda->total_steps - dda->rampup_steps;

#ifdef LOOKAHEAD
    dda->distance = distance;
    dda_find_crossing_speed(prev_dda, dda);
    dda_join_moves(prev_dda, dda);
    dda->n = dda->start_steps;
    if (dda->n == 0)
      dda->c = pgm_read_dword(&c0_P[dda->fast_axis]);
    else
      dda->c = (pgm_read_dword(&c0_P[dda->fast_axis]) *
                int_inv_sqrt(dda->n)) >> 13;
    if (dda->c < dda->c_min)
      dda->c = dda->c_min;
#else
    dda->n = 0;
    dda->c = pgm_read_dword(&c0_P[dda->fast_axis]);
#endif
    dda->c = move_duration / dda->endpoint.F;
    if (dda->c < c_limit)
      dda->c = c_limit;

  } /* ! dda->total_steps == 0 */

  if (DEBUG_DDA && (debug_flags & DEBUG_DDA))
    serial_writestr_P(PSTR("] }\n"));

  // next dda starts where we finish
  memcpy(&startpoint, &dda->endpoint, sizeof(TARGET));
#ifdef LOOKAHEAD
  prev_dda = dda;
#endif
}


void dda_start(DDA *dda) {
  // called from interrupt context: keep it simple!
  if (DEBUG_DDA && (debug_flags & DEBUG_DDA))
    sersendf_P(PSTR("Start: X %lq  Y %lq  Z %lq  F %lu\n"),
               dda->endpoint.axis[X], dda->endpoint.axis[Y],
               dda->endpoint.axis[Z], dda->endpoint.F);

  //sersendf_P(PSTR("dda_started\n"));
  if ( ! dda->nullmove) {
    // get ready to go
    psu_timeout = 0;

    if (dda->endstop_check)
      endstops_on();

    // set direction outputs
    x_direction(dda->x_direction);
    y_direction(dda->y_direction);
    z_direction(dda->z_direction);
    e_direction(dda->e_direction);

    // initialise state variable
    move_state.counter[X] = move_state.counter[Y] = move_state.counter[Z] = \
                            move_state.counter[E] = -(dda->total_steps >> 1);
    memcpy(&move_state.steps[X], &dda->delta[X], sizeof(uint32_t) * 4);
    move_state.endstop_stop = 0;
    move_state.step_no = 0;

    // ensure this dda starts
    dda->live = 1;

    // set timeout for first step
    if (dda->c > 1000) dda->c = 1000;
    timer_set(dda->c, 0);
    //sersendf_P(PSTR("timeout = %lu\n"), dda->c);
  }
  // else just a speed change, keep dda->live = 0

  current_position.F = dda->endpoint.F;
}

void dda_step(DDA *dda) {
  if (move_state.steps[X]) {
    move_state.counter[X] -= dda->delta[X];
    if (move_state.counter[X] < 0) {
      x_step();
      //sersendf_P(PSTR("x"));
      move_state.steps[X]--;
      move_state.counter[X] += dda->total_steps;
    }
  }
  if (move_state.steps[Y]) {
    move_state.counter[Y] -= dda->delta[Y];
    if (move_state.counter[Y] < 0) {
      y_step();
      //sersendf_P(PSTR("y"));
      move_state.steps[Y]--;
      move_state.counter[Y] += dda->total_steps;
    }
  }
  if (move_state.steps[Z]) {
    move_state.counter[Z] -= dda->delta[Z];
    if (move_state.counter[Z] < 0) {
      z_step();
      move_state.steps[Z]--;
      move_state.counter[Z] += dda->total_steps;
    }
  }
  if (move_state.steps[E]) {
    move_state.counter[E] -= dda->delta[E];
    if (move_state.counter[E] < 0) {
      e_step();
      move_state.steps[E]--;
      move_state.counter[E] += dda->total_steps;
    }
  }
  move_state.step_no++;

  // If there are no steps left or an endstop stop happened, we have finished.
  if ((move_state.steps[X] == 0 && move_state.steps[Y] == 0 &&
       move_state.steps[Z] == 0 && move_state.steps[E] == 0) || (move_state.endstop_stop && dda->n <= 0)) {
    dda->live = 0;
    dda->done = 1;
    // sersendf_P(PSTR("end of movement \n"));
#ifdef LOOKAHEAD
    // If look-ahead was using this move, it could have missed our activation:
    // make sure the ids do not match.
    dda->id--;
#endif

    // No need to restart timer here.
    // After having finished, dda_start() will do it.
  }
  else {
    psu_timeout = 0;
    timer_set(dda->c, 0);
  }
  unstep();
}


void dda_clock() {
  DDA *dda;
  static DDA *last_dda = NULL;
  uint8_t endstop_trigger = 0;
  uint32_t move_step_no, move_c;
  int32_t move_n;
  uint8_t recalc_speed;
  uint8_t current_id ;

  dda = queue_current_movement();
  if (dda != last_dda) {
    move_state.debounce_count_x =
      move_state.debounce_count_z =
        move_state.debounce_count_y = 0;
    last_dda = dda;
  }

  if (dda == NULL)
    return;

  if (dda->endstop_check && ! move_state.endstop_stop) {
#ifdef X_MIN_PIN
    if (dda->endstop_check & 0x01) {
      if (x_min() == dda->endstop_stop_cond)
        move_state.debounce_count_x++;
      else
        move_state.debounce_count_x = 0;
      endstop_trigger = move_state.debounce_count_x >= ENDSTOP_STEPS;
    }
#endif

#ifdef Y_MIN_PIN
    if (dda->endstop_check & 0x04) {
      if (y_min() == dda->endstop_stop_cond)
        move_state.debounce_count_y++;
      else
        move_state.debounce_count_y = 0;
      endstop_trigger = move_state.debounce_count_y >= ENDSTOP_STEPS;
    }
#endif


    // If an endstop is definitely triggered, stop the movement.
    if (endstop_trigger) {
      // For always smooth operations, don't halt apruptly,
      // but start deceleration here.
      ATOMIC_START
      move_state.endstop_stop = 1;
      if (move_state.step_no < dda->rampup_steps)  // still accelerating
        dda->total_steps = move_state.step_no * 2;
      else
        // A "-=" would overflow earlier.
        dda->total_steps = dda->total_steps - dda->rampdown_steps +
                           move_state.step_no;
      dda->rampdown_steps = move_state.step_no;
      ATOMIC_END
      // Not atomic, because not used in dda_step().
      dda->rampup_steps = 0; // in case we're still accelerating
      endstops_off();
    }
  } /* ! move_state.endstop_stop */

  ATOMIC_START
  //current_id = dda->id;
  move_step_no = move_state.step_no;
  // All other variables are read-only or unused in dda_step(),
  // so no need for atomic operations.
  ATOMIC_END

  recalc_speed = 0;
  if (move_step_no < dda->rampup_steps) {
#ifdef LOOKAHEAD
    move_n = dda->start_steps + move_step_no;
#else
    move_n = move_step_no;
#endif
    recalc_speed = 1;
  }
  else if (move_step_no >= dda->rampdown_steps) {
#ifdef LOOKAHEAD
    move_n = dda->total_steps - move_step_no + dda->end_steps;
#else
    move_n = dda->total_steps - move_step_no;
#endif
    recalc_speed = 1;
  }
  if (recalc_speed) {
    if (move_n == 0)
      move_c = pgm_read_dword(&c0_P[dda->fast_axis]);
    else
      // Explicit formula: c0 * (sqrt(n + 1) - sqrt(n)),
      // approximation here: c0 * (1 / (2 * sqrt(n))).
      // This >> 13 looks odd, but is verified with the explicit formula.
      move_c = (pgm_read_dword(&c0_P[dda->fast_axis]) *
                int_inv_sqrt(move_n)) >> 13;
    if (move_c < dda->c_min) {
      // We hit max speed not always exactly.
      move_c = dda->c_min;
#ifndef LOOKAHEAD
      dda->rampup_steps = move_step_no;
      dda->rampdown_steps = dda->total_steps - dda->rampup_steps;
#endif
    }
#ifdef LOOKAHEAD
    // Write results.
    ATOMIC_START
    if (current_id == dda->id) {
      dda->c = move_c;
      dda->n = move_n;
    }
    ATOMIC_END
#endif
  }
}

/// update global current_position struct
void update_current_position() {
  DDA *dda = &movebuffer[mb_tail];
  enum axis_e i;

  // Use smaller values to adjust to avoid overflow in later calculations,
  // (STEPS_PER_M_X / 1000) is a bit inaccurate for low STEPS_PER_M numbers.
#ifdef KINEMATICS_CYLINDRICAL
  static const axes_uint32_t PROGMEM steps_per_mm_P = {
    ((STEPS_PER_M_X + 500) / 1000),
    ((STEPS_PER_ROTATION_Y + 500) / 1000),
    ((STEPS_PER_M_Z + 500) / 1000),
    ((STEPS_PER_M_E + 500) / 1000)
  };
#else
  static const axes_uint32_t PROGMEM steps_per_mm_P = {
    ((STEPS_PER_M_X + 500) / 1000),
    ((STEPS_PER_M_Y + 500) / 1000),
    ((STEPS_PER_M_Z + 500) / 1000),
    ((STEPS_PER_M_E + 500) / 1000)
  };
#endif
  if (queue_empty()) {
    for (i = X; i < AXIS_COUNT; i++) {
      current_position.axis[i] = startpoint.axis[i];
    }
  }
  else if (dda->live) {
    for (i = X; i < AXIS_COUNT; i++) {
      current_position.axis[i] = dda->endpoint.axis[i] -
                                 (int32_t)get_direction(dda, i) *
                                 ((move_state.steps[i] * 1000) / pgm_read_dword(&steps_per_mm_P[i]));
    }

    if (dda->endpoint.e_relative)
      current_position.axis[E] =
        (move_state.steps[E] * 1000) / pgm_read_dword(&steps_per_mm_P[E]);

    // current_position.F is updated in dda_start()
  }
}

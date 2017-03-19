
#include	<string.h>
#include	<stdlib.h>
#include	<math.h>
#include	<avr/interrupt.h>
#include    "delay.h"
#include	"dda.h"
#include	"dda_queue.h"
#include "dda_maths.h"
#include	"gcode_parse.h"
#include	"serial.h"			//for debugging
#include	"sermsg.h"			//for debugging
#include	"debug.h"			//for debugging
#include 	"sersendf.h"

#define MIN_SEGMENT 1000

void segmentation(void){
  uint32_t x0 = startpoint.axis[X];
  uint32_t y0 = startpoint.axis[Y];
  uint32_t x1 = next_target.target.axis[X];
  uint32_t y1 = next_target.target.axis[Y];	
  int32_t dx = x1 - x0;
  int32_t dy = y1 - y0;
  
  //length of the move
  uint64_t length = radius_approx(dx,dy);   
  if(length > MIN_SEGMENT){
  uint64_t l = length;
  //find in how many segments we have to divide the move
  uint8_t i = 1;
  while (l > MIN_SEGMENT){
      l = length;
	  l /= i;
	  i++;
}
  TARGET	segment_target;
  segment_target.F = next_target.target.F;
  segment_target.axis[Z] =0;
  // Initialize the extruder axis
  segment_target.axis[E] =0;
  //We start from (x0,y0)
  segment_target.axis[X] = x0;
  segment_target.axis[Y] = y0;
  dx /= i;
  dy /= i;
  for (int k = 0; k < i; k++){
	  segment_target.axis[X] += dx;
	  segment_target.axis[Y] += dy;
	  enqueue(&segment_target);

  }
  }
  next_target.target.axis[Z] = next_target.target.axis[E] = 0;
  enqueue(&next_target.target);
}
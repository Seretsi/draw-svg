#include "viewport.h"

#include "CS248.h"

namespace CS248 {

void ViewportImp::set_viewbox( float x, float y, float span ) {

  // Task 3 (part 2): 
  // Set svg to normalized device coordinate transformation. 
  // Your input arguments are defined as SVG canvans coordinates.

  this->x = x;
  this->y = y;
  this->span = span; 
  Matrix3x3 transform;

  transform(0, 0) = 1 / (2*span); transform(0, 1) = 0; transform(0, 2) = -x/(2*span);
  transform(1, 0) = 0; transform(1, 1) = 1 / (2*span); transform(1, 2) = -y/(2*span);
  transform(2, 0) = 0; transform(2, 1) = 0; transform(2, 2) = 1;
  this->set_canvas_to_norm(transform);
}

void ViewportImp::update_viewbox( float dx, float dy, float scale ) { 
  
  this->x -= dx;
  this->y -= dy;
  this->span *= scale;
  set_viewbox( x, y, span );
}

} // namespace CS248

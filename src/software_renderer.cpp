#include "software_renderer.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>

#include "triangulation.h"

using namespace std;

namespace CS248 {


// Implements SoftwareRenderer //

// fill a sample location with color
void SoftwareRendererImp::fill_sample(int sx, int sy, const Color &color) {
  // Task 2: implement this function

	int ss_width = width * sample_rate;
	int ss_height = height * sample_rate;
	
	// check bounds
	if (sx < 0 || sx >= ss_width) return;
	if (sy < 0 || sy >= ss_height) return;
	if (color.r > 0.0f)
		int i = 0;
	// fill sample - NOT doing alpha blending!
	sample_buffer[4 * (sx + sy * ss_width)] = (uint8_t)(color.r * 255);
	sample_buffer[4 * (sx + sy * ss_width) + 1] = (uint8_t)(color.g * 255);
	sample_buffer[4 * (sx + sy * ss_width) + 2] = (uint8_t)(color.b * 255);
	sample_buffer[4 * (sx + sy * ss_width) + 3] = (uint8_t)(color.a * 255);

}

// fill samples in the entire pixel specified by pixel coordinates
void SoftwareRendererImp::fill_pixel(int x, int y, const Color &color) {

	// Task 2: Re-implement this function

	// check bounds
	if (x < 0 || x >= width) return;
	if (y < 0 || y >= height) return;

	Color pixel_color;
	float inv255 = 1.0f / 255.0f;

	pixel_color.r = pixel_buffer[4 * (x + y * width)] * inv255;
	pixel_color.g = pixel_buffer[4 * (x + y * width) + 1] * inv255;
	pixel_color.b = pixel_buffer[4 * (x + y * width) + 2] * inv255;
	pixel_color.a = pixel_buffer[4 * (x + y * width) + 3] * inv255;

	// todo take sample size into account in mix math
	pixel_color = ref->alpha_blending_helper(pixel_color, color);

	pixel_buffer[4 * (x + y * width)] = (uint8_t)(pixel_color.r * 255);
	pixel_buffer[4 * (x + y * width) + 1] = (uint8_t)(pixel_color.g * 255);
	pixel_buffer[4 * (x + y * width) + 2] = (uint8_t)(pixel_color.b * 255);
	pixel_buffer[4 * (x + y * width) + 3] = (uint8_t)(pixel_color.a * 255);

}

void SoftwareRendererImp::draw_svg( SVG& svg ) {
	// clear sample buffer
	if (!this->sample_buffer.empty())
	{
		this->sample_buffer.assign(this->sample_buffer.size(), 0);
	}
  // set top level transformation
  transformation = canvas_to_screen;

  // canvas outline
  Vector2D a = transform(Vector2D(0, 0)); a.x--; a.y--;
  Vector2D b = transform(Vector2D(svg.width, 0)); b.x++; b.y--;
  Vector2D c = transform(Vector2D(0, svg.height)); c.x--; c.y++;
  Vector2D d = transform(Vector2D(svg.width, svg.height)); d.x++; d.y++;

  svg_bbox_top_left = Vector2D(a.x+1, a.y+1);
  svg_bbox_bottom_right = Vector2D(d.x-1, d.y-1);

  // draw all elements
  for (size_t i = 0; i < svg.elements.size(); ++i) {
    draw_element(svg.elements[i]);
  }

  // draw canvas outline
  rasterize_line(a.x, a.y, b.x, b.y, Color::Black);
  rasterize_line(a.x, a.y, c.x, c.y, Color::Black);
  rasterize_line(d.x, d.y, b.x, b.y, Color::Black);
  rasterize_line(d.x, d.y, c.x, c.y, Color::Black);

  // resolve and send to pixel buffer
  resolve();

}

void SoftwareRendererImp::set_sample_rate( size_t sample_rate ) {

  // Task 2: 
  // You may want to modify this for supersampling support
  this->sample_rate = sample_rate;
  printf("sample rate: %zu\n", sample_rate);

  if (!this->sample_buffer.empty())
  {
	  this->sample_buffer.clear();
  }
  // this->sample_buffer = new unsigned char[4 * this->sample_rate * this->sample_rate * width * height];
  this->sample_buffer.resize(4 * this->sample_rate * this->sample_rate * width * height);
}

void SoftwareRendererImp::set_pixel_buffer( unsigned char* pixel_buffer,
                                             size_t width, size_t height ) {

  // Task 2: 
  // You may want to modify this for supersampling support
  this->pixel_buffer = pixel_buffer;
  this->width = width;
  this->height = height;

  if (!this->sample_buffer.empty())
  {
	  this->sample_buffer.clear();
  }
  this->sample_buffer.resize(4 * this->sample_rate * this->sample_rate * width * height);
}

void SoftwareRendererImp::draw_element( SVGElement* element ) {

	// Task 3 (part 1):
	// Modify this to implement the transformation stack
	 transformation = transformation*element->transform;
	 switch (element->type) {
	case POINT:
		draw_point(static_cast<Point&>(*element));
		break;
	case LINE:
		draw_line(static_cast<Line&>(*element));
		break;
	case POLYLINE:
		draw_polyline(static_cast<Polyline&>(*element));
		break;
	case RECT:
		draw_rect(static_cast<Rect&>(*element));
		break;
	case POLYGON:
		draw_polygon(static_cast<Polygon&>(*element));
		break;
	case ELLIPSE:
		draw_ellipse(static_cast<Ellipse&>(*element));
		break;
	case IMAGE:
		draw_image(static_cast<Image&>(*element));
		break;
	case GROUP:
		draw_group(static_cast<Group&>(*element));
		break;
	default:
		break;
	}
	 transformation = transformation * element->transform.inv();
}


// Primitive Drawing //

void SoftwareRendererImp::draw_point( Point& point ) {

  Vector2D p = transform(point.position);
  rasterize_point( p.x, p.y, point.style.fillColor );

}

void SoftwareRendererImp::draw_line( Line& line ) { 

  Vector2D p0 = transform(line.from);
  Vector2D p1 = transform(line.to);
  rasterize_line( p0.x, p0.y, p1.x, p1.y, line.style.strokeColor );
}

void SoftwareRendererImp::draw_polyline( Polyline& polyline ) {

  Color c = polyline.style.strokeColor;

  if( c.a != 0 ) {
    int nPoints = polyline.points.size();
    for( int i = 0; i < nPoints - 1; i++ ) {
      Vector2D p0 = transform(polyline.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polyline.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_rect( Rect& rect ) {

  Color c;
  
  // draw as two triangles
  float x = rect.position.x;
  float y = rect.position.y;
  float w = rect.dimension.x;
  float h = rect.dimension.y;

  Vector2D p0 = transform(Vector2D(   x   ,   y   ));
  Vector2D p1 = transform(Vector2D( x + w ,   y   ));
  Vector2D p2 = transform(Vector2D(   x   , y + h ));
  Vector2D p3 = transform(Vector2D( x + w , y + h ));
  
  // draw fill
  c = rect.style.fillColor;
  if (c.a != 0 ) {
    rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    rasterize_triangle( p2.x, p2.y, p1.x, p1.y, p3.x, p3.y, c );
  }

  // draw outline
  c = rect.style.strokeColor;
  if( c.a != 0 ) {
    rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    rasterize_line( p1.x, p1.y, p3.x, p3.y, c );
    rasterize_line( p3.x, p3.y, p2.x, p2.y, c );
    rasterize_line( p2.x, p2.y, p0.x, p0.y, c );
  }

}

void SoftwareRendererImp::draw_polygon( Polygon& polygon ) {

  Color c;

  // draw fill
  c = polygon.style.fillColor;
  if( c.a != 0 ) {

    // triangulate
    vector<Vector2D> triangles;
    triangulate( polygon, triangles );

    // draw as triangles
    for (size_t i = 0; i < triangles.size(); i += 3) {
      Vector2D p0 = transform(triangles[i + 0]);
      Vector2D p1 = transform(triangles[i + 1]);
      Vector2D p2 = transform(triangles[i + 2]);
      rasterize_triangle( p0.x, p0.y, p1.x, p1.y, p2.x, p2.y, c );
    }
  }

  // draw outline
  c = polygon.style.strokeColor;
  if( c.a != 0 ) {
    int nPoints = polygon.points.size();
    for( int i = 0; i < nPoints; i++ ) {
      Vector2D p0 = transform(polygon.points[(i+0) % nPoints]);
      Vector2D p1 = transform(polygon.points[(i+1) % nPoints]);
      rasterize_line( p0.x, p0.y, p1.x, p1.y, c );
    }
  }
}

void SoftwareRendererImp::draw_ellipse( Ellipse& ellipse ) {

  // Advanced Task
  // Implement ellipse rasterization

}

void SoftwareRendererImp::draw_image( Image& image ) {

  // Advanced Task
  // Render image element with rotation

  Vector2D p0 = transform(image.position);
  Vector2D p1 = transform(image.position + image.dimension);

  rasterize_image( p0.x, p0.y, p1.x, p1.y, image.tex );
}

void SoftwareRendererImp::draw_group( Group& group ) {

  for ( size_t i = 0; i < group.elements.size(); ++i ) {
    draw_element(group.elements[i]);
  }

}

// Rasterization //

// The input arguments in the rasterization functions 
// below are all defined in screen space coordinates

void SoftwareRendererImp::rasterize_point( float x, float y, Color color ) {

  // fill in the nearest pixel
  int sx = (int)floor(x);
  int sy = (int)floor(y);

  // check bounds
  if (sx < 0 || sx >= width) return;
  if (sy < 0 || sy >= height) return;

  // fill sample - NOT doing alpha blending!
  // TODO: Call fill_pixel here to run alpha blending
  pixel_buffer[4 * (sx + sy * width)] = (uint8_t)(color.r * 255);
  pixel_buffer[4 * (sx + sy * width) + 1] = (uint8_t)(color.g * 255);
  pixel_buffer[4 * (sx + sy * width) + 2] = (uint8_t)(color.b * 255);
  pixel_buffer[4 * (sx + sy * width) + 3] = (uint8_t)(color.a * 255);
}

void swaps(float& x0, float& y0,
	float& x1, float& y1)
{
	auto temp = x0;
	x0 = x1;
	x1 = temp;
	temp = y0;
	y0 = y1;
	y1 = temp;
}

void SoftwareRendererImp::rasterize_line( float x0, float y0,
                                          float x1, float y1,
                                          Color color) {
	// check bounds
	if (x0 < 0 || x0 >= width) return;
	if (y0 < 0 || y0 >= height) return;
	if (x1 < 0 || x1 >= width) return;
	if (y1 < 0 || y1 >= height) return;
	// Task 0: 
  // Implement Bresenham's algorithm (delete the line below and implement your own)
  //ref->rasterize_line_helper(x0, y0, x1, y1, width, height, color, this);
	if (x0 > x1)
		swaps(x0, y0, x1, y1);

	if (y0 <= y1)
	{
		if ((y1 - y0)/(x1 - x0) > 1.0f)
		{
			// mirror the line over y = x
			float temp;
			temp = x0; x0 = y0; y0 = temp;
			temp = x1; x1 = y1; y1 = temp;
			if (x0 > x1) swaps(x0, y0, x1, y1);
			float dx = x1 - x0,
				dy = y1 - y0,
				y = y0,
				eps = 0;

			for (float x = x0; x <= x1; x++) {
				rasterize_point(y, x, color);
				eps += dy;
				if ((eps * 2) >= dx) {
					y++;  eps -= dx;
				}
			}
		}
		else
		{
			//step over x
			float dx = x1 - x0,
				dy = y1 - y0,
				y = y0,
				eps = 0;

			for (float x = x0; x <= x1; x++) {
				rasterize_point(x, y, color);
				eps += dy;
				if ((eps * 2) >= dx) {
					y++;  eps -= dx;
				}
			}
		}
	}
	else
	{
		if ((y1 - y0)/(x1 - x0) < -1.0f)
		{
			// mirror the line over y = x
			float temp;
			temp = x0; x0 = y0; y0 = temp;
			temp = x1; x1 = y1; y1 = temp;
			if (x0 > x1) swaps(x0, y0, x1, y1);
			float dx = x1 - x0,
				dy = y1 - y0,
				y = y0,
				eps = 0;

			for (float x = x0; x <= x1; x++) {
				rasterize_point(y, x, color);
				eps += dy;
				if ((eps * 2) < -dx) {
					y--;  eps += dx;
				}
			}
		}
		else
		{
			//step over x
			float dx = x1 - x0,
				dy = y1 - y0,
				y = y0,
				eps = 0;

			for (float x = x0; x <= x1; x++) {
				rasterize_point(x, y, color);
				eps += dy;
				if ((eps * 2) < -dx) {
					y--;  eps += dx;
				}
			}
		}
	}

  // Advanced Task
  // Drawing Smooth Lines with Line Width
}

void SoftwareRendererImp::rasterize_triangle( float x0, float y0,
                                              float x1, float y1,
                                              float x2, float y2,
                                              Color color ) {
  // Task 1: 
  // Implement triangle rasterization
	//Vector2D aToB(x1 - x0, y1 - y0);
	//Vector2D aToC(x2 - x0, y2 - y0);
	//float p = (float)cross(aToB, aToC);
	//if (p > 0) // epsilon clean up
	//{
	//	float temp;
	//	temp = x1; x1 = x2; x2 = temp;
	//	temp = y1; y1 = y2; y2 = temp;
	//}
	//float dy0 = y1 - y0;
	//float dy1 = y2 - y1;
	//float dy2 = y0 - y2;
	//float dx0 = x1 - x0;
	//float dx1 = x2 - x1;
	//float dx2 = x0 - x2;
	//float c0 = y0 * dx0 - x0 * dy0;
	//float c1 = y1 * dx1 - x1 * dy1;
	//float c2 = y2 * dx2 - x2 * dy2;
	//float minx = min(x0, min(x1, x2));
	//float maxx = max(x0, max(x1, x2));
	//float miny = min(y0, min(y1, y2));
	//float maxy = max(y0, max(y1, y2));
	//for (int x = minx; x <= maxx; x++)
	//	for (int y = miny; y <= maxy; y++)
	//	{
	//		if ((dy0 * (x+0.5f) - dx0 * (y + 0.5f) + c0) < 0.0f) // inside ab
	//			continue;
	//		if ((dy1 * (x + 0.5f) - dx1 * (y + 0.5f) + c2) < 0.0f) // inside ab
	//			continue;
	//		if ((dy2 * (x + 0.5f) - dx2 * (y + 0.5f) + c2) < 0.0f) // inside ab
	//			continue;
	//		rasterize_point(x, y, color);
	//	}
	Vector2D aToB(x1 - x0, y1 - y0);
	Vector2D aToC(x2 - x0, y2 - y0);
	float p = (float)cross(aToB, aToC);
	if (p < 0)
	{
		float temp;
		temp = x1; x1 = x2; x2 = temp;
		temp = y1; y1 = y2; y2 = temp;
	}
	Vector2D ab(x1 - x0, y1 - y0);
	Vector2D bc(x2 - x1, y2 - y1);
	Vector2D ca(x0 - x2, y0 - y2);
	float minx = min(x0, min(x1, x2));
	float maxx = max(x0, max(x1, x2));
	float miny = min(y0, min(y1, y2));
	float maxy = max(y0, max(y1, y2));
	float sample_ctr_offset = 1.0f / (sample_rate * 2.0f);
	float sample_step = 1.0f / sample_rate;

	for (int x = floor(minx); x < maxx; x++)
		for (int y = floor(miny); y < maxy; y++)
		{
			for (int sample_x = 0; sample_x < sample_rate; sample_x++)
			{
				float stride_x = sample_x * sample_step;
				for (int sample_y = 0; sample_y < sample_rate; sample_y++)
				{
					float stride_y = sample_y * sample_step;
					if (((x + stride_x + sample_ctr_offset) - x0) * (y1 - y0) - ((y + stride_y + sample_ctr_offset) - y0) * (x1 - x0) > 0) // inside ab
						continue;
					if (((x + stride_x + sample_ctr_offset) - x1) * (y2 - y1) - ((y + stride_y + sample_ctr_offset) - y1) * (x2 - x1) > 0) // inside ab
						continue;
					if (((x + stride_x + sample_ctr_offset) - x2) * (y0 - y2) - ((y + stride_y + sample_ctr_offset) - y2) * (x0 - x2) > 0) // inside ab
						continue;
					fill_sample(x*sample_rate + sample_x, y*sample_rate + sample_y, color);
				}
			}
		}
	// Advanced Task
  // Implementing Triangle Edge Rules
}

void SoftwareRendererImp::rasterize_image( float x0, float y0,
                                           float x1, float y1,
                                           Texture& tex ) {
  // Task 4: 
  // Implement image rasterization
	for (int x = floor(x0); x < x1; x++)
		for (int y = floor(y0); y < y1; y++)
		{
			float u = (x+0.5 - x0) / (x1+0.5 - x0);
			float v = (y+0.5 - y0) / (y1+0.5 - y0);
			Color c = sampler->sample_nearest(tex, u, v, 0);
			fill_sample(x, y, c);
		}
}

// resolve samples to pixel buffer
void SoftwareRendererImp::resolve( void ) {

  // Task 2: 
  // Implement supersampling
  // You may also need to modify other functions marked with "Task 2".
	// WHAT TO DO: walk through using block grid averaging pixels in ss into image

	for (int x = 0; x < width; x++)
		for (int y = 0; y < height; y++)
		{
			Color c(0.0f, 0.0f, 0.0f, 0.0f);
			int block_sizex = x * sample_rate;
			int block_sizey = y * sample_rate;
			for (int ssx = block_sizex; ssx < block_sizex + sample_rate; ssx++)
				for (int ssy = block_sizey; ssy < block_sizey + sample_rate; ssy++)
				{
					if (4 * (ssx + ssy * width * sample_rate) == 204372)
						int i = 0;
					c.r += ((float)sample_buffer[4 * (ssx + ssy * width * sample_rate)]) / 255.0f;
					c.g += ((float)sample_buffer[4 * (ssx + ssy * width * sample_rate) + 1]) / 255.0f;
					c.b += ((float)sample_buffer[4 * (ssx + ssy * width * sample_rate) + 2]) / 255.0f;
					c.a += ((float)sample_buffer[4 * (ssx + ssy * width * sample_rate) + 3]) / 255.0f;

					if (c.r+c.g+c.b+c.a > 0.0)
						int i = 0;
				}
			c.r /= (float)(sample_rate * sample_rate);
			c.g /= (float)(sample_rate * sample_rate);
			c.b /= (float)(sample_rate * sample_rate);
			c.a /= (float)(sample_rate * sample_rate);
			fill_pixel(x, y, c);
		}
	return;

}

Color SoftwareRendererImp::alpha_blending(Color pixel_color, Color color)
{
  // Task 5
  // Implement alpha compositing
  return pixel_color;
}


} // namespace CS248

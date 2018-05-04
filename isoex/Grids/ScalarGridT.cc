//=============================================================================
//
//  CLASS ScalarGridT - IMPLEMENTATION
//
//=============================================================================

#define ISOEX_SCALARGRIDT_C

//== INCLUDES =================================================================

#include <IsoEx/Grids/ScalarGridT.hh>


//== NAMESPACES ===============================================================

namespace IsoEx {

//== IMPLEMENTATION ==========================================================

float read_float(FILE* _in, bool _swap)
{
  union u3 { float f; unsigned char c[4]; } fc;
  fread((char*)fc.c, 1, 4, _in);
  if (_swap) {
    std::swap(fc.c[0], fc.c[3]);
    std::swap(fc.c[1], fc.c[2]);
  }
  return fc.f;
}

void write_float(float _f, FILE* _out, bool _swap)
{
	union u3 { float f; unsigned char c[4]; } fc;
	fc.f = _f;
	if (_swap) {
		std::swap(fc.c[0], fc.c[3]);
		std::swap(fc.c[1], fc.c[2]);
	}
	fwrite((char*)fc.c, 1, 4, _out);
}


template <class Scalar>
void
ScalarGridT<Scalar>::
sample(const Implicit& _implicit)
{
  for (unsigned int i=0; i<n_points(); ++i)
    values_[i] = _implicit.scalar_distance(point(i));
}


//-----------------------------------------------------------------------------


template <class Scalar>
bool
ScalarGridT<Scalar>::
read(const char* _filename)
{
  bool ok = false;
  FILE* in = fopen(_filename, "rb");
  if (in)
  {
    ok = read(in);
    fclose(in);
  }
  return ok;
}


template <class Scalar>
bool
ScalarGridT<Scalar>::
read(FILE* _in)
{
  // header
  OpenMesh::Vec3f  origin, x_axis, y_axis, z_axis;
  unsigned int     x_res, y_res, z_res;

  
  fscanf(_in, 
	 "%f %f %f\n%f %f %f\n%f %f %f\n%f %f %f\n%d %d %d\n",
	 &origin[0], &origin[1], &origin[2], 
	 &x_axis[0], &x_axis[1], &x_axis[2], 
	 &y_axis[0], &y_axis[1], &y_axis[2], 
	 &z_axis[0], &z_axis[1], &z_axis[2], 
	 &x_res, &y_res, &z_res);

  initialize(origin, x_axis, y_axis, z_axis, x_res, y_res, z_res);
  

  // values
  values_ = Values(n_points(), 0);
  for (unsigned int i=0; i<n_points(); ++i)
    values_[i] = read_float(_in, false);


  return true;
}


//-----------------------------------------------------------------------------


template <class Scalar>
bool
ScalarGridT<Scalar>::
write(const char* _filename)
{
  bool ok = false;
  FILE* out = fopen(_filename, "wb");
  if (out)
  {
    ok = write(out);
    fclose(out);
  }
  return ok;
}


template <class Scalar>
bool
ScalarGridT<Scalar>::
write(FILE* _out)
{
  // header
  fprintf(_out,
	  "%f %f %f\n%f %f %f\n%f %f %f\n%f %f %f\n%d %d %d\n",
	  origin()[0], origin()[1], origin()[2], 
	  x_axis()[0], x_axis()[1], x_axis()[2], 
	  y_axis()[0], y_axis()[1], y_axis()[2], 
	  z_axis()[0], z_axis()[1], z_axis()[2], 
	  x_resolution(), y_resolution(), z_resolution());

  
  // values
  for (unsigned int i=0; i<n_points(); ++i)
	write_float(values_[i], _out, false);


  return true;
}


//=============================================================================
} // namespace IsoEx
//=============================================================================

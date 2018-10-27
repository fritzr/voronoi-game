//------------------------------------------------------------------------------
//  Copyright 2007-2011 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#ifndef _SHPREADER_H_
#define _SHPREADER_H_

#include "polygon.h"
#include "shapelib/shapefil.h"

#include <climits>
#include <string>
#include <vector>
#include <map>

namespace shp
{

template <typename Pt_>
struct PolyData
{
  typedef std::map<unsigned int, std::vector<c_polygon<Pt_> > > polygon_map;
  typedef std::vector<Pt_> point_vector;
  typedef std::map<uint, Pt_> point_map;
};

/* Read a shapefile.
 *
 * Semi-abstract: the default implementation can be used, but will do
 * nothing other than validate the integrity of the file.
 */
template<typename Val_=int>
class ShapeReader
{
private:
  int nShapeType, nShapeFilter, nEntities;
  double adfMinBound[4] = {0.0}, adfMaxBound[4] = { 0.0 };

  SHPHandle shp_;
  DBFHandle dbf_;

  void reset();

protected:
  /* Shape info.  */
  inline int shapeType(void) const { return nShapeType; }
  inline int entities(void) const { return nEntities; }
  inline const double *minBound(void) const { return &adfMinBound[0]; }
  inline const double *maxBound(void) const { return &adfMaxBound[0]; }

  inline SHPHandle shp(void) const { return shp_; }
  inline DBFHandle dbf(void) const { return dbf_; }
  inline size_t rows(void) const { return DBFGetRecordCount(dbf_); }
  inline size_t fields(void) const { return DBFGetFieldCount(dbf_); }

  /* Use this to filter the shapes passed to readShape by type. */
  inline void setFilter(int shapeType) { nShapeFilter = shapeType; }

  /* These are called for you at various points during the shapefile parsing,
   * in the order defined below. If any callbacks return false, the files are
   * closed immediately and no further shapes are read.  */

  /* Called right after the files are opened.
   * All relevant info is available using the getters above.
   * You may also take this opportunity to call setFilter.  */
  virtual bool onOpen(void);

  /* Called for each database field.  */
  virtual bool onField(unsigned int index, DBFFieldType type, const char *name,
      int width, int decimals);

  /* Called for each shape+DBF record in parallel.
   * Use DBFIsAttributeNULL and DBFRead*Attribute to check the DBF fields
   * associated with each row.
   * The shape object is closed for you, do not close it.  */
  virtual bool onRecord(unsigned int row, SHPObject *shp);

  /* Called right before the handles are closed.  */
  virtual void onClose(void);

public:
  ~ShapeReader();
  ShapeReader();

  typedef Val_ value_type;
  virtual value_type value(void) = 0;

  bool read(const std::string& filename);
};

//
// this reader reads ESRI shapefile format and convert
// polygons to poly format
//

template<typename Val_>
class PointReader : public ShapeReader<Val_>
{
public:
  typedef ShapeReader<Val_> super;

protected:
  /* Override.  */
  virtual bool onOpen(void);
  virtual bool onRecord(unsigned int row, SHPObject *shp);

  /* New.  */
  virtual bool onPoint(unsigned int idx, double x, double y) = 0;

  /* Shape info.  */
  using super::shapeType;
  using super::entities;
  using super::minBound;
  using super::maxBound;

  using super::shp;
  using super::dbf;
  using super::rows;
  using super::fields;

public:
  virtual ~PointReader() {}

  PointReader();

  typedef typename ShapeReader<Val_>::value_type value_type;
  virtual value_type value(void) = 0;
};

template<typename Pt_>
class VecPointReader : public PointReader<typename PolyData<Pt_>::point_vector>
{
public:
  typedef typename PolyData<Pt_>::point_vector point_vector;
  typedef PointReader<point_vector> super;
  typedef typename super::value_type value_type;

protected:
  /* Override.  */
  bool onOpen(void) {
    m_points.clear();
    return super::onOpen();
  }

  bool onRecord(unsigned int row, SHPObject *shp)
  {
    if (!super::onRecord(row, shp))
      return false;
    m_points.resize(super::entities());
    return true;
  }

  bool onPoint(unsigned int ptidx, double x, double y)
  {
    m_points.emplace(m_points.begin() + ptidx, Point2d(x, y));
    return true;
  }

public:
  ~VecPointReader() {}
  VecPointReader() : m_points() {}

  /* Implement.  */
  value_type value(void) { return m_points; }

private:
  value_type m_points;
};

template <typename First, typename Second>
std::ostream& operator<<(std::ostream& os, const std::pair<First, Second> &p)
{
  return os << p.first << " => " << p.second;
}

template<typename Pt_>
class IndexedPointReader : public PointReader<typename PolyData<Pt_>::point_map>
{
public:
  typedef typename PolyData<Pt_>::point_map point_map;
  typedef PointReader<point_map> super;
  typedef typename super::value_type value_type;

protected:
  /* Override.  */
  bool onField(unsigned int index, DBFFieldType type, const char *name,
      int width, int decimals);
  bool onPoint(unsigned int ptidx, double x, double y);

public:
  ~IndexedPointReader() {}
  IndexedPointReader() : m_points(), fieldIndex(-1), nextPoint(0u) {}

  /* Implement. */
  value_type value(void) { return m_points; }

private:
  value_type m_points;

  /* Index of the 'pointIndex' field.  */
  int fieldIndex;

  /* If we get a bad pointIndex field, we will just write points
   * incrementally so at least we can get some data. */
  unsigned int nextPoint;
};


template <typename Pt_>
class PolyReader : public ShapeReader<typename PolyData<Pt_>::polygon_map>
{
protected:
  /* Override.  */
  bool onOpen(void);
  bool onField(unsigned int index, DBFFieldType type, const char *name,
      int width, int decimals);
  bool onRecord(unsigned int row, SHPObject *shp);

public:
  typedef typename PolyData<Pt_>::polygon_map polygon_map;
  typedef ShapeReader<polygon_map> super;
  typedef typename super::value_type value_type;

  PolyReader();

  /* Maps pointIndex to a vector of c_polygons.  */
  value_type value(void) { return m_ply_map; }

private:
  /* Index of the 'pointIndex' field.  */
  int fieldIndex;

  /* Index of the 'FTT' field.  */
  int FTTindex;

  /* Last polygon which was inserted, used to add holes.  */
  int lastIndex;

  bool add_ply(int row, c_ply<Pt_>& plys);
  polygon_map m_ply_map;
};

/* Read a shapefile from one of the reader types below.
 * If any errors are encountered, exit the program.  */
template<typename Reader>
typename Reader::value_type
readShapefile(const std::string& path, Reader& r);

template<typename Pt_>
inline typename PolyData<Pt_>::point_vector
readPoints(const std::string& path)
{
  VecPointReader<Pt_> r;
  return readShapefile(path, r);
}

template<typename Pt_>
inline typename PolyData<Pt_>::point_map
readIndexedPoints(const std::string& path)
{
  IndexedPointReader<Pt_> r;
  return readShapefile(path, r);
}

template<typename Pt_>
inline typename PolyData<Pt_>::polygon_map
readPolygons(const std::string& path)
{
  PolyReader<Pt_> r;
  return readShapefile(path, r);
}

} // end namespace shp


// SHP_INSTANTIATE(point_type, [extern])
//
// Declare (with extern) or instantiate (with empty second arg) shape
// reader classes using the given point type.
#define SHP_INSTANTIATE(point_type, E) \
  E template class shp::ShapeReader<\
      typename shp::PolyData<point_type>::point_vector>; \
  E template class shp::ShapeReader<\
      typename shp::PolyData<point_type>::point_map>; \
  E template class shp::ShapeReader<\
      typename shp::PolyData<point_type>::polygon_map>; \
  E template class shp::PointReader<\
      typename shp::PolyData<point_type>::point_vector>; \
  E template class shp::PointReader<\
      typename shp::PolyData<point_type>::point_map>; \
  \
  E template class shp::PolyReader<point_type>; \
  E template class shp::VecPointReader<point_type>; \
  E template class shp::IndexedPointReader<point_type>; \
  \
  E template shp::VecPointReader<point_type>::value_type \
    shp::readShapefile(const std::string& path, VecPointReader<cv::Point2d>& r); \
  E template shp::IndexedPointReader<point_type>::value_type \
    shp::readShapefile(const std::string& path, IndexedPointReader<point_type>& r);\
  E template shp::PolyReader<point_type>::value_type \
    shp::readShapefile(const std::string& path, PolyReader<point_type>& r); \

SHP_INSTANTIATE(cv::Point2d, extern);

#endif //_SHPREADER_H_

//------------------------------------------------------------------------------
//  Copyright 2007-2011 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#ifndef _SHPREADER_H_
#define _SHPREADER_H_

#include "polygon.h"
#include "shapelib/shapefil.h"

#include <climits>
#include <vector>
#include <map>

namespace shp
{

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

  bool read(const string& filename);
};

//
// this reader reads ESRI shapefile format and convert
// polygons to poly format
//

template<typename Val_>
class PointReader : public ShapeReader<Val_>
{
protected:
  /* Override.  */
  virtual bool onOpen(void);
  virtual bool onRecord(unsigned int row, SHPObject *shp);

  /* New.  */
  virtual bool onPoint(unsigned int idx, double x, double y) = 0;

public:
  virtual ~PointReader() {}

  PointReader();

  typedef typename ShapeReader<Val_>::value_type value_type;
  virtual value_type value(void) = 0;
};

typedef std::vector<cv::Point2d> point_vector;

class VecPointReader : public PointReader<point_vector>
{
protected:
  /* Override.  */
  bool onOpen(void) {
    m_points.clear();
    return PointReader::onOpen();
  }

  bool onRecord(unsigned int row, SHPObject *shp)
  {
    if (!PointReader::onRecord(row, shp))
      return false;
    m_points.resize(entities());
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
  typedef typename PointReader<point_vector>::value_type value_type;
  value_type value(void) { return m_points; }

private:
  value_type m_points;
};

template <typename First, typename Second>
std::ostream& operator<<(std::ostream& os, const std::pair<First, Second> &p)
{
  return os << p.first << p.second;
}

typedef std::map<uint, cv::Point2d> point_map;

class IndexedPointReader : public PointReader<point_map>
{
protected:
  /* Override.  */
  bool onField(unsigned int index, DBFFieldType type, const char *name,
      int width, int decimals);
  bool onPoint(unsigned int ptidx, double x, double y);

public:
  ~IndexedPointReader() {}
  IndexedPointReader() : m_points(), fieldIndex(-1), nextPoint(0u) {}

  /* Implement. */
  typedef typename PointReader<point_map>::value_type value_type;
  value_type value(void) { return m_points; }

private:
  value_type m_points;

  /* Index of the 'pointIndex' field.  */
  int fieldIndex;

  /* If we get a bad pointIndex field, we will just write points
   * incrementally so at least we can get some data. */
  unsigned int nextPoint;
};

typedef std::map<unsigned int, std::vector<c_polygon> > polygon_map;

class PolyReader : public ShapeReader<polygon_map>
{
protected:
  /* Override.  */
  bool onOpen(void);
  bool onField(unsigned int index, DBFFieldType type, const char *name,
      int width, int decimals);
  bool onRecord(unsigned int row, SHPObject *shp);

public:
  PolyReader();

  /* Maps pointIndex to a vector of c_polygons.  */
  typedef typename ShapeReader<polygon_map>::value_type value_type;
  value_type value(void) { return m_ply_map; }

private:
  /* Index of the 'pointIndex' field.  */
  int fieldIndex;

  /* Index of the 'FTT' field.  */
  int FTTindex;

  /* Last polygon which was inserted, used to add holes.  */
  int lastIndex;

  bool add_ply(int row, c_ply& plys);
  polygon_map m_ply_map;
};

/* Read a shapefile from one of the reader types below.
 * If any errors are encountered, exit the program.  */
template<typename Reader>
typename Reader::value_type
readShapefile(const string& path, Reader& r);

inline point_vector
readPoints(const string& path)
{
  VecPointReader r;
  return readShapefile(path, r);
}

inline point_map
readIndexedPoints(const string& path)
{
  IndexedPointReader r;
  return readShapefile(path, r);
}

inline polygon_map
readPolygons(const string& path)
{
  PolyReader r;
  return readShapefile(path, r);
}

extern template class ShapeReader<point_vector>;
extern template class PointReader<point_vector>;
extern template class ShapeReader<point_map>;
extern template class PointReader<point_map>;
extern template class ShapeReader<polygon_map>;

extern template VecPointReader::value_type
  readShapefile(const string& path, VecPointReader& r);

extern template IndexedPointReader::value_type
  readShapefile(const string& path, IndexedPointReader& r);

extern template PolyReader::value_type
  readShapefile(const string& path, PolyReader& r);

} // end namespace shp

#endif //_SHPREADER_H_

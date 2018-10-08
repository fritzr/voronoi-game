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

/* Read a shapefile.
 *
 * Semi-abstract: the default implementation can be used, but will do
 * nothing other than validate the integrity of the file.
 */
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
  virtual bool onOpen();

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

  bool read(const string& filename);
};

//
// this reader reads ESRI shapefile format and convert
// polygons to poly format
//

class PointReader : public ShapeReader
{
protected:
  /* Override.  */
  bool onOpen();
  bool onField(unsigned int index, DBFFieldType type, const char *name,
      int width, int decimals);
  bool onRecord(unsigned int row, SHPObject *shp);

public:
  PointReader();

  vector<Point2d>& getPoints() { return m_points; }

private:
  /* Index of the 'pointIndex' field.  */
  int fieldIndex;

  /* If we get a bad pointIndex field, we will just write points
   * incrementally so at least we can get some data. */
  unsigned int nextPoint;

  vector<Point2d> m_points;
};

class PolyReader : public ShapeReader
{
protected:
  /* Override.  */
  bool onOpen();
  bool onField(unsigned int index, DBFFieldType type, const char *name,
      int width, int decimals);
  bool onRecord(unsigned int row, SHPObject *shp);

public:
  typedef typename std::map<unsigned int, vector<c_polygon> > polygon_map;

  /* Maps pointIndex to a vector of c_polygons.  */
  polygon_map& getPolygons(){ return m_ply_map; }

private:
  /* Index of the 'pointIndex' field.  */
  int fieldIndex;

  /* Last polygon which was inserted, used to add holes.  */
  int lastIndex;

  bool add_ply(int row, c_ply& plys, unsigned int index);
  polygon_map m_ply_map;
};

#endif //_SHPREADER_H_

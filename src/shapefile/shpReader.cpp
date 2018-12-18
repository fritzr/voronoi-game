//------------------------------------------------------------------------------
//  Copyright 2007-2011 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#include "shpReader.h"
#include "shapelib/shapefil.h"

#include <cassert>
#include <iostream>
#include <iomanip>
#include <cerrno>

using namespace std;

typedef unsigned int uint;

/* Name of the field which maps polygon/point indexes.  */
static const char POINT_INDEX_FIELD_NAME[] = "pointIndex";

namespace shp
{

  using ::operator<<;

template <typename Val>
ShapeReader<Val>::ShapeReader()
  : nShapeType(SHPT_NULL), nShapeFilter(SHPT_NULL), nEntities(-1),
#ifdef _MSC_VER /* windows */
    adfMinBound({0.0}), adfMaxBound({0.0}),
#else
    adfMinBound{0.0}, adfMaxBound{0.0},
#endif
    shp_(NULL), dbf_(NULL)
{
}

template <typename Val>
ShapeReader<Val>::~ShapeReader()
{
  /* We will just assume our handles are closed properly.  */
}

template <typename Val>
void
ShapeReader<Val>::reset(void)
{
  nShapeType = SHPT_NULL;
  nShapeFilter = SHPT_NULL;
  nEntities = -1;
  memset (&adfMinBound[0], 0, sizeof(adfMinBound));
  memset (&adfMaxBound[0], 0, sizeof(adfMaxBound));

  /* These better be closed already.  */
  shp_ = NULL;
  dbf_ = NULL;
}

template <typename Val>
bool
ShapeReader<Val>::read(const string& filename)
{
  bool ret = true;
  reset();

  shp_ = SHPOpen(filename.c_str(), "rb");
  if (shp_ == NULL)
  {
    cout<<"- Error: Unable to open SHP file "<<filename<<endl;
    return false;
  }

  dbf_ = DBFOpen(filename.c_str(), "rb");
  if (dbf_ == NULL)
  {
    cout<<"- Error: Unable to open DBF file "<<filename<<endl;
    SHPClose(shp_);
    return false;
  }

  // Get info
  SHPGetInfo(shp_, &nEntities, &nShapeType, &adfMinBound[0], &adfMaxBound[0]);

  if (rows() != static_cast<unsigned int>(nEntities))
  {
    cerr << "! Error: Number of database rows (" << rows() << ")"
      " does not match the number of shapes (" << nEntities << ")" << endl;
    goto fail;
  }

  // Open callback
  if (!onOpen())
    goto fail;

  // Visit database fields.
  for (unsigned int fieldidx = 0u; fieldidx < fields(); fieldidx++)
  {
    char fname[12] = { '\0' };
    int width = -1, decs = -1;
    DBFFieldType eType = DBFGetFieldInfo(dbf_, fieldidx, fname, &width, &decs);
    if (!onField(fieldidx, eType, fname, width, decs))
      goto fail;
  }

  // Visit shapes.
  for (unsigned int row = 0u; row < static_cast<uint>(nEntities); row++)
  {
    // Only read non-NULL points that match our filter.
    // If the filter is NULL, "visit" all shapes.
    SHPObject *shape = SHPReadObject(shp_, row);
    if (shape && (nShapeFilter == SHPT_NULL || shape->nSHPType == nShapeFilter))
    {
      if (!onRecord(row, shape))
      {
        SHPDestroyObject(shape);
        goto fail;
      }
    }

    SHPDestroyObject(shape);
  }

  onClose();
  goto success;

fail:
  ret = false;

success:
  DBFClose(dbf_);
  SHPClose(shp_);

  return ret;
}

template <typename Val>
bool
ShapeReader<Val>::onOpen(void)
{
  return true;
}

template <typename Val>
bool
ShapeReader<Val>::onField(unsigned int index, DBFFieldType type,
    const char *name, int width, int decimals)
{
  return true;
}

template <typename Val>
bool
ShapeReader<Val>::onRecord(unsigned int row, SHPObject *shp)
{
  return true;
}

template <typename Val>
void
ShapeReader<Val>::onClose(void)
{
}

/********* PointReader **********/

template<typename Val>
PointReader<Val>::PointReader()
  : ShapeReader<Val>()
{
}

template<typename Val>
bool
PointReader<Val>::onOpen(void)
{
  super::setFilter(SHPT_POINT);
  return true;
}

template<typename Val>
bool
PointReader<Val>::onRecord(unsigned int row, SHPObject *shp)
{
  // Really for SHPT_POINT, there should not be more than 1 vertex,
  // but read all reported vertices anyway.
  for (int vidx = 0; vidx < shp->nVertices; ++vidx)
    if (!onPoint(row, shp->padfX[vidx], shp->padfY[vidx]))
      return false;
  return true;
}

/********* IndexedPointReader *********/

template<typename Pt_>
bool
IndexedPointReader<Pt_>::onField(unsigned int index, DBFFieldType type,
    const char *name, int width, int decimals)
{
  if (0==strncmp(name, POINT_INDEX_FIELD_NAME, sizeof(POINT_INDEX_FIELD_NAME)))
    fieldIndex = index;
  return true;
}

template<typename Pt_>
bool
IndexedPointReader<Pt_>::onPoint(unsigned int row, double x, double y)
{
  if (row == 0)
  {
    if (fieldIndex < 0)
    {
      cerr << "- Error: No DBF field '" << POINT_INDEX_FIELD_NAME << endl;
      return false;
    }
  }

  if (!DBFIsAttributeNULL(super::dbf(), row, fieldIndex))
  {
    int pointIndex = DBFReadIntegerAttribute(super::dbf(), row, fieldIndex);
    if (pointIndex >= 0)
    {
      m_points.emplace(static_cast<unsigned int>(pointIndex), Pt_(x, y));
      nextPoint = pointIndex + 1;
      return true;
    }
    else
    {
      cerr << "- Error: row [" << row << "]: negative pointIndex!" << endl;
      return false;
    }
  }

  /* If there was no known pointIndex, insert counting up from the last known
   * pointIndex (or 0). */
  if (nextPoint < m_points.size())
  {
    cerr << "! Warning: row [" << row << "]: "
      "Automatically inserting at last-known index [" << nextPoint << "]"
      << endl;
    m_points.emplace(nextPoint++, Pt_(x, y));
  }
  else
  {
    cerr << "- Can't find next point index, aborting" << endl;
    return false;
  }

  return true;
}


/********* PolyReader **********/

static const char FTT_FIELD_NAME[] = "FTT";

template<typename Pt_>
PolyReader<Pt_>::PolyReader()
  : fieldIndex(-1), FTTindex(-1), lastIndex(-1), m_ply_map()
{
}

template<typename Pt_>
bool
PolyReader<Pt_>::onOpen(void)
{
  m_ply_map.clear();
  super::setFilter(SHPT_POLYGON);
  return true;
}

template<typename Pt_>
bool
PolyReader<Pt_>::onField(unsigned int index, DBFFieldType type,
    const char *name, int width, int decimals)
{
  if (0==strncmp(name, POINT_INDEX_FIELD_NAME, sizeof(POINT_INDEX_FIELD_NAME)))
    fieldIndex = index;
  else if (0 == strncmp(name, FTT_FIELD_NAME, sizeof(FTT_FIELD_NAME)))
    FTTindex = index;
  return true;
}

template<typename Pt_>
bool
PolyReader<Pt_>::onRecord(unsigned int row, SHPObject *shp)
{
  // Check fields before the first record
  if (row == 0)
  {
    if (fieldIndex < 0)
    {
      cerr << "No DBF field '" << POINT_INDEX_FIELD_NAME << "'!" << endl;
      return false;
    }
    if (FTTindex < 0)
    {
      cerr << "No DBF field '" << FTT_FIELD_NAME << "'!" << endl;
      return false;
    }
  }

  //some format check..
  if (shp->nParts > 0 && shp->panPartStart[0] != 0)
  {
    cerr << "panPartStart[0] = " << shp->panPartStart[0]
      << " not zero as expected" << endl;
    return false;
  }

  //read each part
  int iPart = 1;

  //
  c_ply<Pt_> ply(c_ply<Pt_>::POUT);
  ply.beginPoly();
  ply_extra_info<Pt_>& extra=ply.extra();
  extra.box[0]=shp->dfXMin;
  extra.box[1]=shp->dfXMax;
  extra.box[2]=shp->dfYMin;
  extra.box[3]=shp->dfYMax;
  extra.has_box=true;

  /* Get the index of this polygon's associated center point.  */
  if (!DBFIsAttributeNULL(super::dbf(), row, fieldIndex))
  {
    int pidx = DBFReadIntegerAttribute(super::dbf(), row, fieldIndex);
    if (pidx < 0)
      cerr << "[" << row << "] negative point index!" << endl;
    else
      extra.pointIndex = static_cast<uint>(pidx);
  }

  /* Get the travel time to this ring from the point.  */
  if (!DBFIsAttributeNULL(super::dbf(), row, FTTindex))
  {
    /* This should also be positive.  */
    double ftt = DBFReadDoubleAttribute(super::dbf(), row, FTTindex);
    if (ftt < 0.0)
      cerr << "[" << row << "] negative FTT!" << endl;
    extra.ftt = ftt;
  }

  int j = 0;
  for (; j < shp->nVertices; j++ )
  {
    string pszPartType;

    if( j == 0 && shp->nParts > 0 )
      pszPartType = SHPPartTypeName( shp->panPartType[0] );

    //this defines a new loop...
    if( iPart < shp->nParts && shp->panPartStart[iPart] == j )
    {
      pszPartType = SHPPartTypeName( shp->panPartType[iPart] );
      iPart++;
      ply.endPoly(true); //remove duplicate
      add_ply(j, ply);

      //create a new ply
      ply=c_ply<Pt_>(c_ply<Pt_>::PIN);
      ply.beginPoly();
    }

    //Points may duplicate....
    ply.addVertex(shp->padfX[j],shp->padfY[j], true);
  }

  // use true to set the flag to remove duplicate vertices at the end
  ply.endPoly(true);
  return add_ply(j, ply);
}

template<typename Pt_>
bool
PolyReader<Pt_>::add_ply(int row, c_ply<Pt_>& ply)
{
  if(ply.getType()==c_ply<Pt_>::POUT){
    //create a new polygon
    c_polygon<Pt_> polygon;
    polygon.push_back(ply);
    m_ply_map[ply.extra().pointIndex].push_back(polygon);
    lastIndex = static_cast<int>(ply.extra().pointIndex);
  }

  else
  {
    //c_ply::PIN
    //add to the last polygon
    //make sure that this hole is inside the last polygon...
    if (lastIndex < 0)
    {
      cerr << "! Error: row [" << row << "]: "
        "got a hole as the first layer" << endl;
      return false;
    }

    auto it = m_ply_map.find(static_cast<uint>(lastIndex));
    assert(it != m_ply_map.end());

    c_polygon<Pt_>& polygon = (it->second).back();
    Pt_ pt=ply.findEnclosedPt();
    bool b=polygon.enclosed(pt);
    if(b)
      polygon.push_back(ply);
    else
      cerr<<"! Warning: row [" << row << "]: "
        "hole is ignored due improper nesting"<<endl;
  }
  return true;
}

template<typename Reader>
typename Reader::value_type
readShapefile(const string& path, Reader& r)
{
  bool success = false;
  errno = 0;
  try
  {
    success = r.read(path);
  }
  catch (exception&)
  {
    if (errno)
    {
      cerr << "error opening file '" << path << "': ";
      perror(NULL);
      exit(2);
    }
  }

  if (!success)
  {
    exit(2); // error already reported
  }

#ifdef DEBUG
  auto ref = r.value();
  cerr << "READ shapefile '" << path << "'" << endl;
  unsigned int idx = 0u;
  for (auto it = ref.begin(); it != ref.end(); ++it)
  {
    // may not work for all types
    cerr << "  [" << setw(2) << setfill(' ') << idx++ << "] " << *it << endl;
  }
#endif

  return r.value();
}

} // end namespace shp

SHP_INSTANTIATE(boost::geometry::model::d2::point_xy<double>, );

//------------------------------------------------------------------------------
//  Copyright 2007-2011 by Jyh-Ming Lien and George Mason University
//  See the file "LICENSE" for more information
//------------------------------------------------------------------------------

#include "shpReader.h"
#include "shapelib/shapefil.h"

#include <cassert>
#include <iostream>
using namespace std;

typedef unsigned int uint;

/* Name of the field which maps polygon/point indexes.  */
static const char POINT_INDEX_FIELD_NAME[] = "pointIndex";

ShapeReader::ShapeReader()
  : nShapeType(SHPT_NULL), nShapeFilter(SHPT_NULL), nEntities(-1),
#ifdef _MSC_VER /* windows */
    adfMinBound({0.0}), adfMaxBound({0.0}),
#else
    adfMinBound{0.0}, adfMaxBound{0.0},
#endif
    shp_(NULL), dbf_(NULL)
{
}

ShapeReader::~ShapeReader()
{
  /* We will just assume our handles are closed properly.  */
}

void
ShapeReader::reset(void)
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

bool
ShapeReader::read(const string& filename)
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

bool
ShapeReader::onOpen(void)
{
  return true;
}

bool
ShapeReader::onField(unsigned int index, DBFFieldType type,
    const char *name, int width, int decimals)
{
  return true;
}

bool
ShapeReader::onRecord(unsigned int row, SHPObject *shp)
{
  return true;
}

void
ShapeReader::onClose(void)
{
}

/********* PointReader **********/

PointReader::PointReader()
  : ShapeReader(), fieldIndex(-1), nextPoint(0u), m_points()
{
}

bool
PointReader::onOpen(void)
{
  m_points.clear();
  setFilter(SHPT_POINT);
  return true;
}

bool
PointReader::onField(unsigned int index, DBFFieldType type,
    const char *name, int width, int decimals)
{
  if (0==strncmp(name, POINT_INDEX_FIELD_NAME, sizeof(POINT_INDEX_FIELD_NAME)))
    fieldIndex = index;
  return true;
}

bool
PointReader::onRecord(unsigned int row, SHPObject *shp)
{
  if (row == 0)
  {
    if (fieldIndex < 0)
    {
      cerr << "No DBF field '" << POINT_INDEX_FIELD_NAME << "'!" << endl;
      return false;
    }

    /* Fill with points so we can insert the ordered points at the right
     * location based on the 'pointIndex' field.
     * Practically, these should be in order anyway but there's no telling. */
    m_points.resize(entities());
  }

  bool inserted = false;
  if (!DBFIsAttributeNULL(dbf(), row, fieldIndex))
  {
    unsigned int pointIndex = DBFReadIntegerAttribute(dbf(), row, fieldIndex);
    if (pointIndex < m_points.size())
    {
      m_points[pointIndex] = Point2d(shp->padfX[row], shp->padfY[row]);
      nextPoint = pointIndex + 1;
      inserted = true;
    }
    else
      cerr << "! Warning: row [" << row << "]: pointIndex " << pointIndex
        << " out of range!" << endl;
  }

  if (!inserted)
  {
    if (nextPoint < m_points.size())
    {
      cerr << "! Automatically inserting at last-known index ["
        << nextPoint << "]" << endl;
      m_points[nextPoint++] = Point2d(shp->padfX[row], shp->padfY[row]);
    }
    else
    {
      cerr << "- Can't find next point index, aborting" << endl;
      return false;
    }
  }

  return true;
}

/********* PolyReader **********/

bool
PolyReader::onOpen(void)
{
  m_ply_map.clear();
  setFilter(SHPT_POLYGON);
  return true;
}

bool
PolyReader::onField(unsigned int index, DBFFieldType type,
    const char *name, int width, int decimals)
{
  if (0==strncmp(name, POINT_INDEX_FIELD_NAME, sizeof(POINT_INDEX_FIELD_NAME)))
    fieldIndex = index;
  return true;
}

bool
PolyReader::onRecord(unsigned int row, SHPObject *shp)
{
  if (row == 0)
  {
    if (fieldIndex < 0)
    {
      cerr << "No DBF field '" << POINT_INDEX_FIELD_NAME << "'!" << endl;
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
  c_ply ply(c_ply::POUT);
  ply.beginPoly();
  ply_extra_info& extra=ply.extra();
  extra.box[0]=shp->dfXMin;
  extra.box[1]=shp->dfXMax;
  extra.box[2]=shp->dfYMin;
  extra.box[3]=shp->dfYMax;
  extra.has_box=true;

  int j = 0;
  for (; j < shp->nVertices; j++ )
  {
    string pszPartType;

    /* Get the index of this polygon's associated point.  */
    unsigned int pointIndex = UINT_MAX;
    if (!DBFIsAttributeNULL(dbf(), j, fieldIndex))
      pointIndex = DBFReadIntegerAttribute(dbf(), j, fieldIndex);

    if( j == 0 && shp->nParts > 0 )
      pszPartType = SHPPartTypeName( shp->panPartType[0] );

    //this defines a new loop...
    if( iPart < shp->nParts && shp->panPartStart[iPart] == j )
    {
      pszPartType = SHPPartTypeName( shp->panPartType[iPart] );
      iPart++;
      ply.endPoly(true); //remove duplicate
      add_ply(j, ply, pointIndex);

      //create a new ply
      ply=c_ply(c_ply::PIN);
      ply.beginPoly();
    }

    //Points may duplicate....
    ply.addVertex(shp->padfX[j],shp->padfY[j], true);
  }

  // use true to set the flag to remove duplicate vertices at the end
  ply.endPoly(true);
  return add_ply(j, ply, lastIndex);
}

bool
PolyReader::add_ply(int row, c_ply& ply, unsigned int index)
{
  if(ply.getType()==c_ply::POUT){
    //create a new polygon
    c_polygon polygon;
    polygon.push_back(ply);
    m_ply_map[index].push_back(polygon);
    lastIndex = static_cast<int>(index);
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

    c_polygon& polygon = (it->second).back();
    Point2d pt=ply.findEnclosedPt();
    bool b=polygon.enclosed(pt);
    if(b)
      polygon.push_back(ply);
    else
      cerr<<"! Warning: row [" << row << "]: "
        "hole is ignored due improper nesting"<<endl;
  }
  return true;
}


